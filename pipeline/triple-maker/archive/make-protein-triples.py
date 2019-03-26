#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
import random



__version__="1.0"
__status__ = "Dev"



#####################
def random_string(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))



###########################
def load_json_db(data_dir):

	db_obj = {}
        data_dir = "/data/projects/glygen/generated/output/"

        file_list1 = glob.glob(data_dir + "human_protein_*.csv")
        file_list2 = glob.glob(data_dir + "mouse_protein_*.csv")
        file_list3 = glob.glob(data_dir + "*_proteoform_glycosylation_sites_unicarbkb_glytoucan.csv")
                        
        exclude_list = [ 
                "human_protein_isoform2ensembl.csv"
        ]

        fasta_files = glob.glob(data_dir + "*_protein_all.fasta")
        seq_hash = {}
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = record.id.split("|")[1]
                desc = record.description
                seq_hash[seq_id] = str(record.seq.upper())
       

        tmpl_hash = {}
        with open("outdir/mapping.txt", 'r') as FR:
            data_frame = csv.reader(FR, delimiter='|', quotechar='"')
            rowCount = 0
            for row in data_frame:
                rowCount += 1
                tmpl_hash[row[0]] = row[1]



        FW = open("outdir/protein-fields.txt", "w")
        for inFile in file_list1 + file_list2:
                fileName = os.path.basename(inFile)
                if fileName in exclude_list:
                    continue
                prefix = "_".join(fileName.split(".")[0].split("_")[2:])
                with open(inFile, 'r') as FR:
			dataFrame = csv.reader(FR, delimiter=',', quotechar='"')
			rowCount = 0
			field_list = []
			for row in dataFrame:
                                rowCount += 1
				if rowCount == 1:
					field_list = row
				        FW.write(">%s %s\n%s\n\n" % (fileName,prefix,"\n".join(row)))
                                else:
					main_id_index = field_list.index("uniprotkb_acc_canonical")
                                        id1 = row[main_id_index]
					if id1 not in db_obj:
						db_obj[id1] = {}
					if prefix not in db_obj[id1]:
						db_obj[id1][prefix] = []
					row_obj = {}
                                        for j in xrange(1,len(field_list)):
                                            id2 = field_list[j]
					    row_obj[id2] = [] if row[j].strip() == "" else row[j].split("|")
                                            if "" in row_obj[id2]:
                                                row_obj[id2].remove("")
                                        db_obj[id1][prefix].append(row_obj)
        FW.close()

        for inFile in file_list3:
                fileName = os.path.basename(inFile)
                if fileName in exclude_list:
                    continue
                prefix = "_".join(fileName.split(".")[0].split("_")[2:])
                with open(inFile, 'r') as FR:
                        dataFrame = csv.reader(FR, delimiter=',', quotechar='"')
                        rowCount = 0
                        field_list = []
                        for row in dataFrame:
                                rowCount += 1
                                if rowCount == 1:
                                        field_list = row
                                else:
                                        main_id_index = field_list.index("uniprotkb_acc_canonical")
                                        id1 = row[main_id_index]
                                        if id1 not in db_obj:
                                                db_obj[id1] = {}
                                        if prefix not in db_obj[id1]:
                                                db_obj[id1][prefix] = []
                                        row_obj = {}
                                        for j in xrange(0,len(field_list)):
                                                id2 = field_list[j]
                                                row_obj[id2] = [] if row[j].strip() == "" else row[j].split("|")
                                                if "" in row_obj[id2]:
                                                    row_obj[id2].remove("") 
                                        db_obj[id1][prefix].append(row_obj)

        return db_obj, seq_hash, tmpl_hash




#######################################
def main():


        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]


        taxid2name = {9606:"Homo sapiens", 10090:"Mus musculus"}
        data_dir = "/data/projects/glygen/generated/output/"
        db_obj, seq_hash, tmpl_hash = load_json_db(data_dir)

        out_obj_list = []
	for id1 in db_obj:
            if "idmapping" in  db_obj[id1]:

                #Extract accession
                uniprotkb_ac = id1.split("-")[0]

                #Extract from "idmapping"
                gene = []
                gene_symbol = ""
                if len(db_obj[id1]["idmapping"][0]["gene_name"]) > 0:
                    gene_symbol = db_obj[id1]["idmapping"][0]["gene_name"][0]
                hgnc_url = ""
                if "uniprot_all_idmapping" in db_obj[id1]:
                    if len(db_obj[id1]["uniprot_all_idmapping"][0]["hgnc_id"]) > 0:
                        hgnc_url = root_obj["glygen_api_gene_url"] 
                        hgnc_url += db_obj[id1]["uniprot_all_idmapping"][0]["hgnc_id"][0]
                if gene_symbol != "":
                    gene.append({"name":gene_symbol
                        ,"url":hgnc_url
                    })

                #Extract from "information"
                species = []
                keywords = []
                protein_mass, protein_length, tax_id,  uniprotkb_id = -1.0, -1, -1, ""
                if "information" in db_obj[id1]:
                    if len(db_obj[id1]["information"][0]["keywords_label"]) > 0:
                        keywords = db_obj[id1]["information"][0]["keywords_label"]
                    if len(db_obj[id1]["information"][0]["protein_mass"]) > 0:
                        protein_mass = float(db_obj[id1]["information"][0]["protein_mass"][0])
                    if len(db_obj[id1]["information"][0]["protein_length"]) > 0:
                        protein_length = int(db_obj[id1]["information"][0]["protein_length"][0])
                    if len(db_obj[id1]["information"][0]["tax_id"]) > 0:
                        tax_id = int(db_obj[id1]["information"][0]["tax_id"][0])
                        species += [
                            {
                                "name":taxid2name[tax_id]
                                ,"taxid":tax_id
                                ,"evidence":[
                                    {
                                        "database":"UniProtKB"
                                        ,"id":uniprotkb_ac
                                        ,"url":"http://uniprot.org"
                                    }
                                ]
                            }
                        ]
                    if len(db_obj[id1]["information"][0]["uniprotkb_id"]) > 0:
                        uniprotkb_id = db_obj[id1]["information"][0]["uniprotkb_id"][0]
                #Extract from "recnames"
                recommended_name = {}
                if "recnames" in db_obj[id1]:
                    if len(db_obj[id1]["recnames"][0]["recommended_name_full"]) > 0:
                        recommended_name["full"] = db_obj[id1]["recnames"][0]["recommended_name_full"][0]
                    if len(db_obj[id1]["recnames"][0]["recommended_name_short"]) > 0:
                        recommended_name["short"] = db_obj[id1]["recnames"][0]["recommended_name_short"][0]

                #Extract from "alternativename"
                alternative_names = []
                if "alternativename" in db_obj[id1]:
                    for o in db_obj[id1]["alternativename"]:
                        full_name, short_name = "", ""
                        if len(o["alternative_name_full"]) > 0:
                            full_name = o["alternative_name_full"][0]
                        if len(o["alternative_name_short"]) > 0:
                            short_name = o["alternative_name_short"][0]
                        alternative_names += [{"full":full_name, "short":short_name}]


                
                #Extract from "glycosylation"
                glycosylation = []
                if "glycosylation_sites_unicarbkb_glytoucan" in db_obj[id1]:
                    for obj in db_obj[id1]["glycosylation_sites_unicarbkb_glytoucan"]:
                        glycan_id = obj["glytoucan_acc"][0]
                        r_val = random_string(4, "LVFMCAGPSTYWQNHIEDKR")
                        glycan_type = random_string(4, "LVFMCAGPSTYWQNHIEDKR")
                        o = {
                                "glycan_id":glycan_id
                                ,"glycan_type":glycan_type
                                ,"position":int(obj["canonical_glycosylation_site"][0])
                                ,"residue":obj["amino_acid"][0]
                                ,"evidence":[]
                        }
                        for ev in obj["evidence"]:
                            r_val = random_string(4, "LVFMCAGPSTYWQNHIEDKR")
                            o["evidence"].append({"database":r_val, "id":r_val, "url":r_val})
                        glycosylation.append(o)

                #Extract from "pathway"
                pathway = []
                if "pathway_accessions" in db_obj[id1]:
                    for obj in db_obj[id1]["pathway_accessions"]:
                        for p_id in obj["kegg_id"]:
                            pathway.append({"id":p_id, "name":"", "resource":"KEGG"})
                        for p_id in obj["reactome_id"]:
                            pathway.append({"id":p_id, "name":"", "resource":"Reactome"})

                #Extract from "crossref"
                crossref = []
                r_val = random_string(4, "LVFMCAGPSTYWQNHIEDKR")
                crossref.append({"database":r_val, "id":r_val, "url":r_val})

                #Extract sequence for canonical
                sequence = {
                    "sequence":seq_hash[id1]
                    ,"length":len(seq_hash[id1])
                }

                #Extract isoforms
                isoforms = []
                for obj in db_obj[id1]["idmapping"]:
                    for isoform in obj["reviewed_isoforms"] + obj["unreviewed_isoforms"]:
                        r_val = random_string(4, "LVFMCAGPSTYWQNHIEDKR")
                        isoform_seq = seq_hash[isoform] if isoform in seq_hash else r_val
                        isoforms.append({
                            "id":isoform
                            ,"url":root_obj["uniprot_isoform_url"] + isoform
                            ,"sequence":{
                                "sequence":isoform_seq
                                ,"length":len(isoform_seq)
                            }
                        })
                        t_key,s_val = "isoform_type", isoform,
                        triple = tmpl_hash[t_key] % (s_val)
                        print triple

                        t_key,s_val,o_val = "isoform_status", isoform, str(isoform in obj["reviewed_isoforms"]).lower()
                        triple = tmpl_hash[t_key] % (s_val, o_val)
                        print triple

                        t_key,s_val,o_val = "isoform_canonical", isoform, isoform == id1
                        triple = tmpl_hash[t_key] % (s_val, o_val)
                        print triple

                        t_key,s_val,o_val = "isoform_mass", isoform, 0.0
                        triple = tmpl_hash[t_key] % (s_val, o_val)
                        print triple



                publication = []
                r_val = random_string(4, "LVFMCAGPSTYWQNHIEDKR")
                publication.append({"pmid":r_val, "title":r_val, "journal":r_val,"url":r_val})



                out_obj = {
                    "accession":id1
                    ,"id":uniprotkb_id
                    ,"mass":protein_mass
                    ,"recommendedname":recommended_name
                    ,"alternativenames":alternative_names
                    ,"species":species
                    ,"gene":gene
                    ,"crossref":crossref
                    ,"sequence":sequence
                    ,"isoforms":isoforms
                    ,"glycosylation":glycosylation
                    ,"keywords":keywords
                    ,"publication":publication
                    ,"pathway":pathway
                }
                
                out_obj_list.append(out_obj)








if __name__ == '__main__':
        main()



