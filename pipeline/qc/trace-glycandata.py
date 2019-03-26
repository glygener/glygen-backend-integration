import os,sys
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

from Bio import SeqIO
import pymongo
from pymongo import MongoClient
import glob


__version__="1.0"
__status__ = "Dev"

######################
def load_dataframe(data_frame, sheet_name, in_file, separator):

    field_list = []
    data_frame[sheet_name] = {}
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='"')
        row_count = 0
        for row in csv_grid:
            row_count += 1
            if row_count == 1:
                field_list = row
            else:
                row_obj = {}
                for j in xrange(1,len(field_list)):
                    field_name  = field_list[j]
                    row_obj[field_name] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                main_id = row[0].strip()
                if main_id not in data_frame[sheet_name]:
                    data_frame[sheet_name][main_id] = []
                data_frame[sheet_name][main_id].append(row_obj)
    return field_list



##############################
def get_recnames():

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    
    ac2canon = {}

    query_string = """
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?ac ?isoformuri
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:sequence ?isoformuri .
            ?isoformuri gly:canonical 'true' ^^<http://www.w3.org/2001/XMLSchema#boolean> .
        }
        """
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        canon = result["isoformuri"]["value"].split("/")[-1]
        ac2canon[ac] = canon


    sparql.setQuery("""
        PREFIX rdf:   <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX gly: <http://purl.glygen.org/>
        PREFIX up: <http://purl.uniprot.org/core/>
        SELECT ?ac ?fullname
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:recommendedName ?nameuri .
            ?nameuri up:fullName ?fullname . 
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
        }
    """)
    
    recname_dict = {}
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        full_name = result["fullname"]["value"].strip()
        recname_dict[ac] = full_name

        #if ac in ac2canon:
	#    canon = ac2canon[ac]
	#    recname_dict[canon] = full_name

    return recname_dict



###############################
def main():

		
	usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

        file_list = [ 
                "downloads/glytoucan/gtc_v3/mass.txt"
                ,"downloads/glytoucan/gtc_v3/allmotif.txt"
                ,"downloads/glytoucan/gtc_v3/uckb.txt"
                ,"downloads/glytoucan/gtc_v3/pbch.txt"
                ,"downloads/glytoucan/gtc_v3/gdb.txt"
                ,"downloads/glytoucan/gtc_v3/enzyme/gtc2enz.txt"
                ,"downloads/glytoucan/gtc_v3/taxa.txt"
        ]
        recname_dict = get_recnames()


        canon2taxid = {}
        with open("reviewed/mouse_protein_idmapping.csv", "r") as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            for row in csv_grid:
                canon2taxid[row[0]] = 9606
                ac = row[0].split("-")[0]
                canon2taxid[ac] = 9606


	src_grid = {}
        for in_file in file_list:
            with open(in_file, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter='\t', quotechar='"')
                row_count = 0
                for row in csv_grid:
                    row_count += 1
                    if row_count == 1:
                        field_list = row
                    else:
                        row_obj = {}
                        for j in xrange(1,len(field_list)):
                            field_name  = field_list[j]
                            row_obj[field_name] = [] if row[j].strip() == "" else [row[j].replace("\"","")]
                            #row_obj[field_name] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                        main_id = row[0].strip()
                        if main_id not in src_grid:
                            src_grid[main_id] = []
                        src_grid[main_id].append(row_obj)

        seq_dict = {"glycoct":{}, "iupac":{}, "wurcs":{}}
        for db in ["iupac", "wurcs", "glycoct"]:
            for f in glob.glob("downloads/glytoucan/gtc_v3/%s/*.txt" % (db)):
                ac = os.path.basename(f)[0:-4]
                seq = ""
                with open(f, "r") as FR:
                    for line in FR:
                        seq += line.strip() + " "
                seq_dict[db][ac] = seq.strip()


        extra_files = [
            "downloads/uckb/human_glycosylation.csv"
            ,"downloads/uckb/mouse_glycosylation.csv"
        ]
        seen_gsite = {}
        for in_file in extra_files:
            with open(in_file, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
                for row in csv_grid:
                    ac = row[0].split("/")[-1]
                    pos = row[1].split("^^")[0]
                    uckb_id = row[3]
                    combo_id = "%s,%s,%s" %(ac,pos,uckb_id)
                    seen_gsite[combo_id] = True


        val_dict = {}
        for main_id in src_grid:
            if main_id not in val_dict:
                val_dict[main_id] = {}
            for row in src_grid[main_id]:
                for field in row:
                    if field not in val_dict[main_id]:
                        val_dict[main_id][field] = []
                    val_dict[main_id][field] += row[field]



        client = MongoClient('mongodb://localhost:27017')
	db = client["glygen"]
        
        cond_list = []
        #cond_list.append({"species.taxid": {'$eq': int(tax_id)}})
        #cond_list.append({"glycosylation.evidence.database": {'$eq':"UniCarbKB"}})
        query_obj = {} if cond_list == [] else { "$and": cond_list }
        for doc in db["c_glycan"].find(query_obj):
            ac = doc["glytoucan_ac"]
            
            if True: #check mass and Monosaccharides count
                mass_one, mass_two = doc["mass"], round(float(src_grid[ac][0]["DerivatizedMass"][0]), 2)
                flag = mass_one == mass_two
                print "mass|%s|%s" % (ac,flag)

                val_one, val_two = doc["number_monosaccharides"], round(float(src_grid[ac][0]["Monosaccharides"][0]), 2)
                flag = val_one == val_two
                print "mono_count|%s|%s" % (ac,flag)


            if True: #Check motif
                val_one, val_two = [], []
                for obj in doc["motifs"]:
                    val_one += obj.values()
                for obj in src_grid[ac]:
                    if "MotifAccession" in obj:
                        val_two += [obj["MotifAccession"][0], obj["Label"][0]]
                val_one = sorted(set(val_one))
                val_two = sorted(set(val_two))
                flag = val_one == val_two
                print "motif|%s|%s|%s|%s" % (ac,flag,val_one, val_two)

            if True: #Check sequences
                for field in ["iupac", "wurcs", "glycoct"]:
                    val_one = seq_dict[field][ac] if ac in seq_dict[field] else ""
                    val_two = doc[field]
                    flag = val_one == val_two
                    print "%s|%s|%s" % (field,ac,flag)
                    
            if True: #Check cross ref
                val_one, val_two = [], []
                for obj in doc["crossref"]:
                    val_two += [obj["id"]]
                for obj in src_grid[ac]:
                    for f in ["PubChem","GlycomeDBID","UniCarbKBID"]:
                        if f in obj:
                            val_one += [obj[f][0]]
                val_one = sorted(set(val_one))
                val_two = sorted(set(val_two))
                flag = val_one == val_two
                print "crossref|%s|%s|%s|%s" % (ac,flag,val_one, val_two)


            if True: #check enzyme (for mouse only since human enzymes are derived through orthology)
                tmp_val = [[], [], []]
                for obj in doc["enzyme"]:
                    if obj["uniprot_canonical_ac"] in canon2taxid:
                        tmp_val[0].append(obj["uniprot_canonical_ac"].split("-")[0])
                        tmp_val[1].append(obj["gene"])
                        tmp_val[2].append(obj["protein_name"])            
                field_list = ["Enzyme:UniProtAccession", "Enzyme:GeneName", "Enzyme:Description"]

                val_dict[ac]["Enzyme:Description"] = []
                if "Enzyme:UniProtAccession" in val_dict[ac]:
                    for enzyme_ac in val_dict[ac]["Enzyme:UniProtAccession"]:
                        if enzyme_ac in recname_dict:
                            val_dict[ac]["Enzyme:Description"].append(recname_dict[enzyme_ac])
                for i in xrange(0, len(field_list)):
                    field = field_list[i]
                    if field not in val_dict[ac]:
                        val_dict[ac][field] = []
                    if "" in val_dict[ac][field]:
                        val_dict[ac][field].remove("")
                    if "" in tmp_val[i]:
                        tmp_val[i].remove("")
                    flag = sorted(set(val_dict[ac][field])) == sorted(set(tmp_val[i]))
                    print "%s|%s|%s" % (field, ac,flag)


            if True: #check classficiation
                val_one, val_two = [], []
                for obj in doc["classification"]:
                    val_two.append(obj["type"]["name"] + "|" + obj["subtype"]["name"]);
                class_list = [
                    "N-Glycan complex"
                    ,"N-Glycan hybrid"
                    ,"N-Glycan high mannose"
                    ,"O-Glycan core 1"
                    ,"O-Glycan core 2"
                    ,"O-Glycan core 3"
                    ,"O-Glycan core 4"
                    ,"O-Glycan core 5"
                    ,"O-Glycan core 6"
                    ,"O-Glycan core 7"
                ]
                if "Label" in val_dict[ac]:
                    for val in val_dict[ac]["Label"]:
                        if val in class_list:
                            gtype = val.split(" ")[0]
                            gsubtype = " ".join(val.split(" ")[1:])
                            val_one.append(gtype + "|" + gsubtype);
                flag = sorted(set(val_one)) == sorted(set(val_two))
                print "classification|%s|%s|%s|%s" % (ac,flag, val_one, val_two)



            if True: #check species
                tmp_val = [[], []]
                for obj in doc["species"]:
                    tmp_val[0].append(str(obj["taxid"]))
                if "UniProtTaxID" not in val_dict[ac]:
                    val_dict[ac]["UniProtTaxID"] = []
                for val in val_dict[ac]["UniProtTaxID"]:
                    if val in ["9606", "10090"]:
                        tmp_val[1].append(val)
                flag = sorted(set(tmp_val[0])) == sorted(set(tmp_val[1]))
                print "species|%s|%s" % (ac,flag)
                #Check glycoprotein (since mapping glytoucan_ac to uckb_id is already validated
                #it suffices to make sure uniprot_ac+pos+uckb_id is seen in csv downloaded from uckb
                
            if True: #check glycosylation
                tv_list = []
                for obj in doc["glycoprotein"]:
                    t_list = []
                    for o in obj["evidence"]:
                        uckb_id = o["id"]
                        combo_id = "%s,%s,%s" % (obj["uniprot_canonical_ac"].split("-")[0], 
                                obj["position"],uckb_id) 
                        tv_list.append(combo_id in seen_gsite)
                if len(tv_list) > 0:
                    print "glycoprotein", ac, sorted(set(tv_list))[0]


if __name__ == '__main__':
	main()


