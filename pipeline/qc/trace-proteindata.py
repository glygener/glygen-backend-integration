import os,sys
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


######################
def load_dataframe(data_frame, sheet_name, in_file, separator):

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
    return





##########################
def load_ac2xreflist(ep_data, tax_id):


    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    query_string = """
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>      
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?ac ?xrefuri  ?dbname
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac rdfs:seeAlso ?xrefuri .
            ?xrefuri up:database ?dbname . 
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
        }
        LIMIT 100000000
        """

    query_string = query_string % (tax_id)
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    seen_pmid = {}
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1].strip()
        db_id = result["xrefuri"]["value"].split("/")[-1].strip()
        db_name = result["dbname"]["value"].split("/")[-1].strip()
        if ac not in ep_data["ac2xreflist"]:
            ep_data["ac2xreflist"][ac] = []
        ep_data["ac2xreflist"][ac].append({"dbid":db_id, "dbname":db_name})

        print ac, db_id, db_name
    sys.exit()


    return



##########################
def load_canon2citation(ep_data, tax_id):

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    query_string = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?ac ?citeuri ?journalname ?journaltitle
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:citation ?citeuri .
            ?citeuri up:name ?journalname .
            ?citeuri up:title ?journaltitle .
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
        }
        """

    query_string = query_string % (tax_id)
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    seen_pmid = {}
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1].strip()
        canon = ep_data["ac2canon"][ac] if ac in ep_data["ac2canon"] else ""
        if canon != "":
            pm_id = result["citeuri"]["value"].split("/")[-1]
            combo_id = "%s %s" % (canon, pm_id)
            if combo_id not in seen_pmid:
                seen_pmid[combo_id] = True
                journal_name = result["journalname"]["value"].strip().replace("\"", "")
                journal_title = result["journaltitle"]["value"].strip().replace("\"", "")
                if canon not in ep_data["canon2citelist"]:
                    ep_data["canon2citelist"][canon] = []
                o = {"pmid":pm_id, "journalname":journal_name, "journaltitle":journal_title}
                ep_data["canon2citelist"][canon].append(o)

    return




##########################
def load_canon2genelist(ep_data, tax_id):

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    query_string = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?ac ?genename
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:encodedBy ?gene_uri .
            ?gene_uri skos:prefLabel ?genename .
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
        }
        """

    query_string = query_string % (tax_id)
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        canon = ep_data["ac2canon"][ac] if ac in ep_data["ac2canon"] else ""
        if canon != "":
            gene_name = result["genename"]["value"].strip()
            if canon not in ep_data["canon2genelist"]:
                ep_data["canon2genelist"][canon] = []
            ep_data["canon2genelist"][canon].append(gene_name)

    return




##########################
def load_isoform2mass(ep_data, tax_id):

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    query_string = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?isoform ?mass
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?isoform up:mass ?mass .
            ?isoform rdf:type <http://purl.uniprot.org/core/Simple_Sequence> .
        }
        """
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        isoform = result["isoform"]["value"].split("/")[-1]
        mass = result["mass"]["value"].split("/")[-1]
        ep_data["isoform2mass"][isoform] = mass
    return


#######################
def load_ac2fullname(ep_data, tax_id):

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
        
    query_string = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX gly: <http://purl.glygen.org/>
        PREFIX up: <http://purl.uniprot.org/core/>
        SELECT ?ac ?fullname
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:recommendedName ?nameuri .
            ?nameuri up:fullName ?fullname . 
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
        }
    """
    
    query_string = query_string % (tax_id)
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        ep_data["ac2fullname"][ac] = result["fullname"]["value"].strip()
    return


#######################
def load_ac2shortname(ep_data, tax_id):

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
        
    query_string = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX gly: <http://purl.glygen.org/>
        PREFIX up: <http://purl.uniprot.org/core/>
        SELECT ?ac ?shortname
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:recommendedName ?nameuri .
            ?nameuri up:shortName ?shortname .
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
        }
    """
    
    query_string = query_string % (tax_id)
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        ep_data["ac2shortname"][ac] = result["shortname"]["value"].strip()
    return


##########################
def load_isoform2seq(ep_data, tax_id):



    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    query_string = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?isoformuri ?isoformseq ?iscanon
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?isoformuri gly:canonical ?iscanon .
            ?isoformuri rdf:value ?isoformseq .
        }
        """
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        isoform = result["isoformuri"]["value"].split("/")[-1]
        seq = result["isoformseq"]["value"].strip()
        is_canon = result["iscanon"]["value"]
        ep_data["isoform2seq"][isoform] = seq

    return


##########################
def load_isoform2locusinfo(ep_data, tax_id):

    #?trsid gly:reverseStrand ?strand .
    #?trsid gly:transcriptRange ?rangeuri .

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    query_string = """
        PREFIX faldo: <http://biohackathon.org/resource/faldo#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?isoformuri ?trsid ?pepid ?chrid ?startpos ?endpos
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?isoformuri gly:ensTranscript ?trsid .
            ?trsid up:translatedTo ?pepid . 
            ?trsid gly:chromosome ?chrid .
            ?trsid gly:transcriptRange ?rangeuri .
            ?rangeuri faldo:begin ?beginuri .
            ?rangeuri faldo:end ?enduri .
            ?beginuri faldo:position ?startpos .
            ?enduri faldo:position ?endpos .
        }
        """

    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        isoform = result["isoformuri"]["value"].split("/")[-1]
        if isoform not in ep_data["isoform2locusinfo"]:
            ep_data["isoform2locusinfo"][isoform] = {
                    "chrid":result["chrid"]["value"].strip()
                    ,"trsid":result["trsid"]["value"].split("/")[-1] 
                    ,"pepid":result["pepid"]["value"].split("/")[-1]
                    ,"strand":"x"
                    ,"startpos":result["startpos"]["value"]
                    ,"endpos":result["endpos"]["value"]
            }


    return



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--configfile",action="store",dest="configfile",help="config JSON file")
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

    (options,args) = parser.parse_args()
    for file in ([options.species, options.configfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    config_obj = json.loads(open(options.configfile, "r").read())
    species = options.species
    tax_id = "9606" if species == "human" else "10090"


    client = MongoClient('mongodb://localhost:27017')
    db = client["glygen"]

    ep_data = {
        "ac2canon":{}, 
        "ac2fullname":{}, 
        "ac2shortname":{},
        "canon2genelist":{},
        "canon2citelist":{},
        "isoform2seq":{},
        "isoform2mass":{},
        "isoform2locusinfo":{},
        "ac2xreflist":{}
    }

    
    #Over write ac2canon from mongodb for now
    query_obj = {}
    for doc in db["c_protein"].find(query_obj):
        canon = doc["uniprot_canonical_ac"]
        ac = canon.split("-")[0]
        ep_data["ac2canon"][ac] = canon 

    #load_isoform2mass(ep_data, tax_id)
    #load_isoform2locusinfo(ep_data, tax_id)
    #load_canon2citation(ep_data, tax_id)
    #load_ac2fullname(ep_data, tax_id)
    #load_ac2shortname(ep_data, tax_id)
    #load_canon2genelist(ep_data, tax_id)
    
    #load_isoform2seq(ep_data, tax_id)
    #load_ac2xreflist(ep_data, tax_id)


    data_frame = {}
    if species == "human":
        #in_file = config_obj["pathinfo"]["downloads"] + "/biomuta/biomuta.csv"
        #load_dataframe(data_frame, "mutation", in_file, ",")
        
        #in_file = config_obj["pathinfo"]["downloads"] + "/bioxpress/human_protein_expression_disease.2.csv"
        #load_dataframe(data_frame, "expression_disease", in_file, ",")
        
        in_file = config_obj["pathinfo"]["downloads"] + "/bioxpress/human_protein_expression_normal.2.csv"
        load_dataframe(data_frame, "expression_normal", in_file, ",")



    client = MongoClient('mongodb://localhost:27017')
    db = client["glygen"]
    #query_obj = {}
    query_obj = {"species.taxid": {'$eq': int(tax_id)}}

    for doc in db["c_protein"].find(query_obj):
        canon = doc["uniprot_canonical_ac"]
        ac = canon.split("-")[0]
        if False: #check fullname, shortname
            val_one, val_two = "", ""
            if ac in ep_data["ac2fullname"]:
                val_one =  ep_data["ac2fullname"][ac]
            if "full" in doc["recommendedname"]:
                val_two = doc["recommendedname"]["full"]
            flag = val_one ==  val_two
            print "fullname|%s|%s|%s|%s|%s" % (flag, ac, canon, val_one, val_two)
        
            val_one, val_two = "", ""
            if ac in ep_data["ac2shortname"]:
                val_one = ep_data["ac2shortname"][ac]
            if "short" in doc["recommendedname"]:
                val_two = doc["recommendedname"]["short"]
            flag = val_one ==  val_two
            print "shortname|%s|%s|%s|%s" % (flag, canon, val_one, val_two)
                        
        if False: #check mass
            val_one = int(ep_data["isoform2mass"][canon]) if canon in ep_data["isoform2mass"] else ""
            val_two = int(doc["mass"]["chemical_mass"])
            flag =  val_one == val_two
            print "mass|%s|%s|%s|%s" % (flag, canon, val_one, val_two)

        if False: #check gene names
            val_one, val_two = [], []
            for obj in doc["gene"]:
                val_two.append(obj["name"])
            if canon in ep_data["canon2genelist"]:
                val_one = ep_data["canon2genelist"][canon]
            flag = val_one == val_two
            print "gene_name|%s|%s|%s|%s" % (flag, canon, val_one, val_two)

        if False: #Check publications 
            val_one, val_two = [], []
            if canon in ep_data["canon2citelist"]:
                for obj in ep_data["canon2citelist"][canon]:
                    val_one += obj.values()
            for obj in doc["publication"]:
                obj.pop("url")
                val_two += obj.values()
            val_one = sorted(set(val_one))
            val_two = sorted(set(val_two))
            for i in xrange(0, len(val_one)):
                val_one[i] = str(val_one[i]).lower()
            for i in xrange(0, len(val_two)):
                val_two[i] = str(val_two[i]).lower()
            flag = val_one == val_two
            print "publications|%s|%s|%s|%s" % (flag, canon, len(val_one), len(val_two))
            #if canon == "Q58FF3-1":
            #    print val_one
            #    print val_two

        if False: #Check cross references
            val_one_a, val_two_a = [], []
            val_one_b, val_two_b = [], []
            if ac in ep_data["ac2xreflist"]:
                for obj in ep_data["ac2xreflist"][ac]:
                    val_one_a.append(obj["dbid"])
                    val_one_b.append(obj["dbname"])
            for obj in doc["crossref"]:
                val_two_a.append(obj["id"])
                val_two_b.append(obj["database"])

            val_one_a = sorted(set(val_one_a))
            val_two_a = sorted(set(val_two_a))
            val_one_b = sorted(set(val_one_b))
            val_two_b = sorted(set(val_two_b))

            flag = set(val_one_a) & set(val_two_a) == set(val_two_a)
            print "crossref_dbid|%s|%s|%s|%s" % (flag, canon, val_one_a, val_two_a)
            flag = set(val_one_b) & set(val_two_b) == set(val_two_b)
            print "crossref_dbname|%s|%s|%s|%s" % (flag, canon, val_one_b, val_two_b)


        if False: #check mutation
            val_one, val_two = [], []
            if canon in data_frame["mutation"]:
                for obj in data_frame["mutation"][canon]:
                    val_one += [obj["aa_pos"][0], obj["ref_aa"][0],obj["alt_aa"][0],obj["do_id"][0]]
            for obj in doc["mutation"]:
                val_two += [str(obj["start_pos"]),str(obj["sequence_org"]),str(obj["sequence_mut"]),str(obj["disease"]["doid"])]
            val_one = sorted(set(val_one))
            val_two = sorted(set(val_two))
            flag = set(val_one) & set(val_two) == set(val_two)
            print "mutation|%s|%s|%s|%s" % (flag, canon, val_one, val_two)

        if False: #check expression_disease
            val_one, val_two = [], []
            if ac in data_frame["expression_disease"]:
                for obj in data_frame["expression_disease"][ac]:
                    val_one += [obj["direction"][0], obj["significance"][0],obj["do_id"][0]]
            for obj in doc["expression_disease"]:
                val_two += [str(obj["trend"]),str(obj["significant"]),str(obj["disease"]["doid"])]
            val_one = sorted(set(val_one))
            val_two = sorted(set(val_two))
            flag = set(val_one) & set(val_two) == set(val_two)
            print "expression_disease|%s|%s|%s|%s" % (flag, canon, val_one, val_two)


        if True: #check expression_normal
            val_one, val_two = [], []
            if ac in data_frame["expression_normal"]:
                for obj in data_frame["expression_normal"][ac]:
                    call = "yes" if obj["expressionCall"][0] == "present" else "no"
                    val_one += [obj["uberonAnatomyId"][0].replace("UBERON:", ""), call]
            for obj in doc["expression_tissue"]:
                val_two += [str(obj["tissue"]["uberon"]),str(obj["present"])]
            val_one = sorted(set(val_one))
            val_two = sorted(set(val_two))
            flag = set(val_one) & set(val_two) == set(val_two)
            print "expression_normal|%s|%s|%s|%s" % (flag, canon, val_one, val_two)



        if False: #Check canon sequence
            val_one, val_two = "", ""
            if canon in ep_data["isoform2seq"]:
                val_one = ep_data["isoform2seq"][canon]
            if "sequence" in doc["sequence"]:
                val_two = doc["sequence"]["sequence"]
            flag = val_one == val_two
            print "canon_seq|%s|%s|%s|%s" % (flag, canon, val_one, val_two)


        if False: #Check isoform information
            for obj in doc["isoforms"]:
                isoform = obj["isoform_ac"]
                val_one, val_two = "", ""
                if isoform in ep_data["isoform2seq"]:
                    val_one = ep_data["isoform2seq"][isoform]
                if "sequence" in obj["sequence"]:
                    val_two = obj["sequence"]["sequence"]
                if val_one != "": #only check if endpoint gave sequence
                    flag = val_one == val_two
                    print "isoform_seq|%s|%s|%s|%s" % (flag, isoform, val_one, val_two)

                val_one, val_two = "", ""
                if isoform in ep_data["isoform2locusinfo"]:
                    val_one = ep_data["isoform2locusinfo"][isoform]["chrid"]
                val_two = obj["locus"]["chromosome"]
                flag = val_one == val_two
                print "chromosome|%s|%s|%s|%s" % (flag, isoform, val_one, val_two)

                val_one, val_two = "", ""
                if isoform in ep_data["isoform2locusinfo"]:
                    val_one = ep_data["isoform2locusinfo"][isoform]["trsid"]
                if len(obj["locus"]["evidence"]) > 0:
                    val_two =  obj["locus"]["evidence"][0]["id"]
                flag = val_one == val_two
                print "trsid|%s|%s|%s|%s" % (flag, isoform, val_one, val_two)

                val_one, val_two = "", ""
                if isoform in ep_data["isoform2locusinfo"]:
                    val_one = ep_data["isoform2locusinfo"][isoform]["pepid"]
                if len(obj["locus"]["evidence"]) > 0:
                    val_two =  obj["locus"]["evidence"][1]["id"]
                flag = val_one == val_two
                print "pepid|%s|%s|%s|%s" % (flag, isoform, val_one, val_two)


                val_one, val_two = 0, 0
                if isoform in ep_data["isoform2locusinfo"]:
                    val_one = ep_data["isoform2locusinfo"][isoform]["startpos"]
                if "start_pos" in obj["locus"]:
                    val_two = obj["locus"]["start_pos"]
                flag = int(val_one) == int(val_two)
                print "startpos|%s|%s|%s|%s" % (flag, isoform, val_one, val_two) 
            
                val_one, val_two = 0, 0
                if isoform in ep_data["isoform2locusinfo"]:
                    val_one = ep_data["isoform2locusinfo"][isoform]["endpos"]
                if "end_pos" in obj["locus"]:
                    val_two = obj["locus"]["end_pos"]
                flag = int(val_one) == int(val_two)
                print "endpos|%s|%s|%s|%s" % (flag, isoform, val_one, val_two)  
             
            




if __name__ == '__main__':
	main()
