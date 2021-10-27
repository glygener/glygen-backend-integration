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
import commands



sys.path.append('../../glytools/')
import libgly




def load_abstracts(abstract_dict):

    file_list = glob.glob("downloads/ncbi/medline/pmid.*.txt")
    total = len(file_list)
    f_count = 0
    for in_file in file_list:
        pmid = in_file.split(".")[-2]
        FR = open(in_file, "r")
        line_list = []
        tag = ""
        for line in FR:
            tag = line[0:5] if line[0:4].strip() != "" else tag
            if tag == "AB  -":
                line_list.append(line[5:].strip())
        if line_list != []:
            combo_id = "pubmed."+pmid
            abstract_dict[combo_id] = " ".join(line_list).strip()
        FR.close()
        f_count += 1


    return



def extract_name(obj_list, name_type, resource):
   
    if obj_list == []:
        return ""
    
    name_list_dict = {"recommended":[], "synonym":[]}
    for obj in obj_list:
        if obj["resource"] == resource:
            name_list_dict[obj["type"]].append(obj["name"])


    
    if name_type == "recommended" and name_list_dict["recommended"] == []:
        name_list_dict["recommended"] += name_list_dict["synonym"]


    if name_type == "all":
        return "; ".join(name_list_dict["recommended"] + name_list_dict["synonym"])
    else:
        return "; ".join(name_list_dict[name_type])





def get_protein_record(obj, record_dict):

    canon = obj["uniprot_canonical_ac"]
    mass = -1

            

    protein_name, protein_names_uniprotkb, protein_names_refseq = "", "", ""
    if "protein_names" in obj:
        protein_name = extract_name(obj["protein_names"], "recommended", "UniProtKB")
        protein_names_uniprotkb = extract_name(obj["protein_names"], "all", "UniProtKB")
        protein_names_refseq = extract_name(obj["protein_names"], "all", "RefSeq")

    refseq_ac, refseq_name = "", ""
    if "refseq" in obj:
        refseq_ac = obj["refseq"]["ac"] if "ac" in obj["refseq"] else ""
        refseq_name = obj["refseq"]["name"] if "name" in obj["refseq"] else ""

    gene_name, gene_names_uniprotkb, gene_names_refseq = "", "", ""
    if "gene_names" in obj:
        gene_name = extract_name(obj["gene_names"], "recommended", "UniProtKB")
        gene_names_uniprotkb = extract_name(obj["gene_names"], "all", "UniProtKB")
        gene_names_refseq = extract_name(obj["gene_names"], "all", "RefSeq")

    
    org_obj_list = []
    if "species" in obj:
        for o in obj["species"]:
            o["common_name"] = species_obj[str(o["taxid"])]["common_name"]
            org_obj_list.append(o)



            
    publication_count = 0
    if "publication" in obj:
        for pub_obj in obj["publication"]:
            for o in pub_obj["reference"]:
                record_id = ""
                if type(o["type"]) is list:
                    o["type"] = o["type"][0]
                record_id = "%s.%s" % (o["type"].lower(), o["id"])
                if record_id not in record_dict:
                    record_dict[record_id] = pub_obj
                    record_dict[record_id]["referenced_proteins"] = [canon]
                    record_dict[record_id]["referenced_glycans"] = []
                    if "species" not in record_dict[record_id]:
                        record_dict[record_id]["species"] = []
                    record_dict[record_id]["species"] += org_obj_list
                else:
                    record_dict[record_id]["referenced_proteins"].append(canon)


    for sec in site_sec_list:
        if sec in obj:
            for sec_obj in obj[sec]:
                for ev_obj in sec_obj["evidence"]:
                    record_id = ""
                    if ev_obj["database"] != "":
                        record_id = "%s.%s" % (ev_obj["database"].lower(), ev_obj["id"])
                    if record_id != "" and record_id in record_dict:
                        if sec not in record_dict[record_id]:
                            record_dict[record_id][sec] = []
                        sec_obj["uniprot_canonical_ac"] = canon
                        record_dict[record_id][sec].append(sec_obj)



    return 





def get_glycan_record(doc, record_dict, seen):

    glytoucan_ac = doc["glytoucan_ac"] 

    org_obj_list = []
    if "species" in doc:
        for o in doc["species"]:
            o["common_name"] = species_obj[str(o["taxid"])]["common_name"]
            org_obj_list.append(o)



    sec = "publication"
    if sec in doc:
        for sec_obj in doc[sec]:
            for o in sec_obj["reference"]:
                record_id = ""
                if type(o["type"]) is list:
                    o["type"] = o["type"][0]
                record_id = "%s.%s" % (o["type"].lower(),o["id"])
                if record_id not in record_dict:
                    record_dict[record_id] = sec_obj
                    record_dict[record_id]["referenced_proteins"] = []
                    record_dict[record_id]["referenced_glycans"] = [glytoucan_ac]
                    if "species" not in record_dict[record_id]:
                        record_dict[record_id]["species"] = [] 
                    record_dict[record_id]["species"] += org_obj_list
                else:
                    record_dict[record_id]["referenced_glycans"].append(glytoucan_ac)

   
    sec = "expression" 
    if sec in doc:
        for sec_obj in doc[sec]:
            if sec_obj["tissue"] == {} and sec_obj["cell_line"] == {} and sec_obj["abundance"] == "":
                continue
            expr_obj = {
                "uniprot_canonical_ac":sec_obj["uniprot_canonical_ac"],
                "residue":sec_obj["residue"],
                "start_pos":sec_obj["start_pos"],
                "end_pos":sec_obj["end_pos"],
                "glytoucan_ac": glytoucan_ac,
                "abundance":sec_obj["abundance"],
                "tissue":sec_obj["tissue"],
                "cell_line":sec_obj["cell_line"]
            }
            expr_obj_str = json.dumps(expr_obj)

            for ev_obj in sec_obj["evidence"]:
                record_id = ""
                if "database" in ev_obj:
                    record_id = "%s.%s" % (ev_obj["database"].lower(),ev_obj["id"]) 
                if record_id != "" and record_id in record_dict:
                    if record_id not in seen:
                        seen[record_id] = {}
                    
                    if "glcan_expression" not in record_dict[record_id]:
                        record_dict[record_id]["glycan_expression"] = []
                    if expr_obj != {} and expr_obj_str not in seen[record_id]:
                        record_dict[record_id]["glycan_expression"].append(expr_obj)
                        seen[record_id][expr_obj_str] = True




    return




def main():

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global main_dict
    global site_sec_list

    config_file = "../../conf/config-1.1.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    site_sec_list = config_obj["sitesections"]

    abstract_dict = {}
    load_abstracts(abstract_dict)

    species_obj = {}
    in_file = path_obj["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)


    record_dict = {}
    seen = {}
    for record_type in ["protein", "glycan"]:
        #file_list = glob.glob("jsondb/%sdb/G05049IC.json" % (record_type))
        #file_list = glob.glob("jsondb/%sdb/P14210*.json" % (record_type))
        #file_list = glob.glob("jsondb/%sdb/Q8K4H4*.json" % (record_type))
        file_list = glob.glob("jsondb/%sdb/*.json" % (record_type))

        record_count = 0
        for json_file in file_list:
            doc = json.loads(open(json_file,"r").read())
            if record_type == "protein":
                get_protein_record(doc, record_dict)
            elif record_type == "glycan":
                get_glycan_record(doc, record_dict,seen)
    
    record_count = 0
    for record_id in record_dict:
        record_obj = record_dict[record_id]
        record_obj["reference"] = record_obj["reference"][0] 
        record_obj["record_id"] = record_id
        ab = ""
        if record_id.lower().find("pubmed") != -1 and record_id in abstract_dict:
            ab = abstract_dict[record_id]

        record_obj["abstract"] = ab
        if "evidence" in record_obj:
            record_obj.pop("evidence")
        for sec in site_sec_list:
            if sec not in record_obj:
                record_obj[sec] = []
        
        out_file = "jsondb/publicationdb/%s.json" % (record_id.lower().replace("/", "_"))
        out_str = json.dumps(record_obj, indent=4)
        if len(out_str) > 16000000:
            out_file = "jsondb/jumbodb/%s.json" % (record_id.lower().replace("/", "_"))

        with open(out_file, "w") as FW:
            FW.write("%s\n" % (out_str))
        record_count += 1
    print "make-publicationdb: ... final created: %s publication records" % (record_count)
                                                                                 
 




if __name__ == '__main__':
    main()

