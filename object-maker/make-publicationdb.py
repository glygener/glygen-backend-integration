#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
import subprocess
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq



import libgly
import csvutil
import batchutil
import section_stats



def get_uncited_pubs():


    is_cited = {}
    file_list = glob.glob("reviewed/*_citations*.csv")
    for in_file in file_list:
        if in_file.find(".stat.csv") != -1:
            continue
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            if len(row) != len(f_list):
                print (in_file)
                print (f_list)
                print (row)
                exit()
            xref_key, xref_id = row[f_list.index("xref_key")], row[f_list.index("xref_id")]
            if xref_key.split("_")[-1] not in ["pubmed","doi"]:
                continue
            is_cited[xref_id] = True


    file_list = glob.glob("generated/misc/*_sectioninfo.json")
    seen_ptrn = {}
    for in_file in file_list:
        doc = json.loads(open(in_file, "r").read())
        for sec in doc:
            for ptrn in doc[sec]["sheetlist"]:
                seen_ptrn[ptrn] = True

    row_list = []
    seen = {}
    for ptrn in seen_ptrn:
        file_list = glob.glob("reviewed/*%s*" % (ptrn))
        for in_file in file_list:
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            f_list = data_frame["fields"]
            if "xref_key" not in f_list:
                continue
            for row in data_frame["data"]:
                xref_key, xref_id = row[f_list.index("xref_key")], row[f_list.index("xref_id")]
                if xref_key.split("_")[-1] not in ["pubmed","doi"]:
                    continue
                c = "%s,%s" % (xref_id, in_file)
                if c in seen:
                    continue
                seen[c] = True
                if xref_id not in is_cited:
                    row_list.append([xref_id, in_file])

    return row_list




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
    tax_name = ""
    if "species" in obj:
        for o in obj["species"]:
            o["common_name"] = species_obj[str(o["taxid"])]["common_name"]
            org_obj_list.append(o)
        if org_obj_list != []:
            tax_name = obj["species"][0]["common_name"]

    r_protein = {"uniprot_canonical_ac":canon, "protein_name":protein_name, "gene_name":gene_name, "tax_name":tax_name }
                    
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
                    record_dict[record_id]["referenced_glycans"] = {}
                    record_dict[record_id]["referenced_proteins"] = {canon:r_protein}
                    if "species" not in record_dict[record_id]:
                        record_dict[record_id]["species"] = []
                    record_dict[record_id]["species"] += org_obj_list
                else:
                    record_dict[record_id]["referenced_proteins"][canon] = r_protein


    for sec in site_sec_list:
        if sec in ["biomarkers"]:
            continue
        if sec in obj:
            for sec_obj in obj[sec]:
                if "evidence"  not in sec_obj:
                    continue
                for ev_obj in sec_obj["evidence"]:
                    record_id = ""
                    if ev_obj["database"].lower() in ["doi", "pubmed"]:
                        record_id = "%s.%s" % (ev_obj["database"].lower(), ev_obj["id"])
                    if record_id != "" and record_id in record_dict:
                        if sec not in record_dict[record_id]:
                            record_dict[record_id][sec] = []
                        sec_obj["uniprot_canonical_ac"] = canon
                        record_dict[record_id][sec].append(sec_obj)

    #print (list(record_dict.keys()))
    sec = "biomarkers"
    load_biomarker_records(sec, obj[sec], record_dict, "uniprot_canonical_ac", canon, "protein")  
    return 



def load_biomarker_records(sec, sec_obj_list, record_dict, main_id_field, main_id, record_type):    

    for sec_obj in sec_obj_list:
        sec_obj[main_id_field] = main_id
        for ev_obj in sec_obj["evidence"]:
            record_id = ""
            if ev_obj["database"].lower() in ["doi", "pubmed"]:
                record_id = "%s.%s" % (ev_obj["database"].lower(), ev_obj["id"])
            if record_id != "" and record_id in record_dict:
                if sec not in record_dict[record_id]:
                    record_dict[record_id][sec] = {}
                o = {"evidence":[]}
                for f in ["assessed_biomarker_entity", "biomarker_id", "biomarker", "condition"]:
                    o[f] = sec_obj[f]
                for oo in sec_obj["evidence"]:
                    o["evidence"].append(oo)
                combo_id = "%s" % (o["biomarker_id"])
                if combo_id not in record_dict[record_id][sec]:
                    record_dict[record_id][sec][combo_id] = o
                else:
                    for oo in sec_obj["evidence"]:
                        record_dict[record_id][sec][combo_id]["evidence"].append(oo)
        


    return



def get_glycan_record(doc, record_dict, seen):




    glytoucan_ac = doc["glytoucan_ac"] 

    org_obj_list = []
    if "species" in doc:
        for o in doc["species"]:
            o["common_name"] = species_obj[str(o["taxid"])]["common_name"]
            org_obj_list.append(o)
            #for oo in o["evidence"]:
            #    if oo["database"] in ["PubMed", "DOI"]:
            #        record_id = "%s.%s" % (oo["database"].lower(),oo["id"])
            #        print (record_id)



    sec = "publication"
    if sec in doc:
        for sec_obj in doc[sec]:
            for o in sec_obj["reference"]:
                record_id = ""
                if type(o["type"]) is list:
                    o["type"] = o["type"][0]
                record_id = "%s.%s" % (o["type"].lower(),o["id"])
                #print (record_id)
                if record_id not in record_dict:
                    record_dict[record_id] = sec_obj
                    record_dict[record_id]["referenced_proteins"] = {}
                    record_dict[record_id]["referenced_glycans"] = {glytoucan_ac:True}
                    if "species" not in record_dict[record_id]:
                        record_dict[record_id]["species"] = [] 
                    record_dict[record_id]["species"] += org_obj_list
                else:
                    record_dict[record_id]["referenced_glycans"][glytoucan_ac] = True


    sec = "biomarkers"
    load_biomarker_records(sec, doc[sec], record_dict, "glytoucan_ac", glytoucan_ac, "glycan")

 
    sec = "expression" 
    if sec in doc:
        for sec_obj in doc[sec]:
            if sec_obj["tissue"] == {} and sec_obj["cell_line"] == {} and sec_obj["abundance"] == "":
                continue
            expr_obj = {
                "glytoucan_ac": glytoucan_ac,
                "abundance":sec_obj["abundance"],
                "tissue":sec_obj["tissue"],
                "cell_line":sec_obj["cell_line"]
            }
            for f in ["uniprot_canonical_ac", "residue", "start_pos", "end_pos"]:
                if f in sec_obj:
                    expr_obj[f] = sec_obj[f]

            expr_obj_str = json.dumps(expr_obj)
            for ev_obj in sec_obj["evidence"]:
                record_id = ""
                if "database" in ev_obj:
                    record_id = "%s.%s" % (ev_obj["database"].lower(),ev_obj["id"]) 
                if record_id != "" and record_id in record_dict:
                    if record_id not in seen:
                        seen[record_id] = {}
                    if "glycan_expression" not in record_dict[record_id]:
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

    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    site_sec_list = config_obj["sitesections"]

    url_map = {
        "pubmed":"https://pubmed.ncbi.nlm.nih.gov/%s",
        "doi":"https://doi.org/%s"
    }

    #DEBUG = True
    DEBUG = False



    log_file = "logs/make-publicationdb.log"
    msg = "make-publicationdb: stareted logging"
    csvutil.write_log_msg(log_file, msg, "w")

    #glob_pattern_list = ["G23863VK*.json"]
    glob_pattern_list = ["P14210*.json", "G17689DH*.json"]

    abstract_dict = {}
    if DEBUG == False:
        #Clear existing jumbo objects of this record type
        cmd = "rm -f jsondb/jumbodb/publicationdb/*"
        x = subprocess.getoutput(cmd)
        cmd = "rm -f jsondb/batchdb/publication.*"
        x = subprocess.getoutput(cmd)
        glob_pattern_list = ["*.json"]
        #row_list = get_uncited_pubs()
        #if row_list != []:
        #    print ("The following IDs are not in *_citations_* datasets")
        #    print (json.dumps(row_list, indent=4))
        #    exit()
        msg = "make-publicationdb: loading abstracts"
        csvutil.write_log_msg(log_file, msg, "a")
        load_abstracts(abstract_dict)
        

    species_obj = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)



    record_dict = {}
    seen = {}
    for record_type in ["protein", "glycan"]:
        file_list = []
        for glob_pattern in glob_pattern_list:
            file_list += glob.glob("jsondb/%sdb/%s" % (record_type, glob_pattern))

        record_count = 0
        for json_file in file_list:
            main_id = json_file.split("/")[-1].replace(".json", "")
            if record_count > 0 and record_count%1000 == 0:
                msg = "make-publicationdb: processed %s %s records" % (record_count, record_type)
                csvutil.write_log_msg(log_file, msg, "a")

            jumbo_file = "jsondb/jumbodb/%sdb/%s.json" % (record_type, main_id)
            if os.path.isfile(jumbo_file):
                json_file = jumbo_file
            doc = json.loads(open(json_file,"r").read())
            if record_type == "protein":
                get_protein_record(doc, record_dict)
                record_count += 1
            elif record_type == "glycan":
                get_glycan_record(doc, record_dict,seen)
                record_count += 1
        msg = "make-publicationdb: processed %s %s records" % (record_count, record_type)
        csvutil.write_log_msg(log_file, msg, "a")

     
    # For publication objects, references to PubMed and DOI
    # should point to external pages instead of GlyGen publication page
    for record_id in record_dict:
        for obj in record_dict[record_id]["reference"]:
            prefix = record_id.split(".")[0]
            xref_id = ".".join(record_id.split(".")[1:])
            if prefix in url_map:
                obj["url"] = url_map[prefix] % (xref_id)
        for sec in site_sec_list:
            if sec in ["biomarkers"]:
                continue
            if sec in record_dict[record_id]:
                for obj in record_dict[record_id][sec]:
                    for o in obj["evidence"]:
                        if o["database"] in ["DOI", "PubMed"]:
                            prefix = o["database"].lower()
                            o["url"] = url_map[prefix] % (o["id"])


    record_count = 0
    for record_id in record_dict:
        record_obj = record_dict[record_id]
        tmp_list = []
        for ac in record_obj["referenced_glycans"]:
            tmp_list.append({"glytoucan_ac":ac})
        record_obj["referenced_glycans"] = tmp_list
        tmp_list = []
        for canon in record_obj["referenced_proteins"]:
            tmp_list.append(record_obj["referenced_proteins"][canon])
        record_obj["referenced_proteins"] = tmp_list
        record_obj["reference"] = record_obj["reference"][0] 
        record_obj["record_id"] = record_id
        ab = ""
        if record_id.lower().find("pubmed") != -1 and record_id in abstract_dict:
            ab = abstract_dict[record_id]

        record_obj["abstract"] = ab
        if "evidence" in record_obj:
            record_obj.pop("evidence")
        for sec in site_sec_list:
            if sec in ["biomarkers"]:
                continue
            if sec not in record_obj:
                record_obj[sec] = []
        if "biomarkers" in record_obj:
            tmp_list = []
            for biomarker_id in record_obj["biomarkers"]:
                o = record_obj["biomarkers"][biomarker_id]
                #seen_evdn = {}
                #evdn_list = []
                #for oo in o["evidence"]:
                #    oo_str = json.dumps(oo)
                #    if oo_str not in seen_evdn:
                #        evdn_list.append(oo)
                #        seen_evdn[oo_str] = True
                #o["evidence"] = evdn_list
                if "evidence" in o:
                    o.pop("evidence")
                tmp_list.append(o)
            record_obj["biomarkers"] = tmp_list

        section_stats.get_sec_stats(record_obj, "publication")

        out_file = "jsondb/publicationdb/%s.json" % (record_id.lower().replace("/", "_"))
        out_str = json.dumps(record_obj, indent=4)
        if len(out_str) > 16000000:
            out_file = "jsondb/jumbodb/publicationdb/%s.json" % (record_id.lower().replace("/", "_"))
        if DEBUG:
            print ("created %s" % (out_file)) 
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (out_str))
        record_count += 1
    msg = "make-publicationdb: ... final created: %s publication records" % (record_count)
    csvutil.write_log_msg(log_file, msg, "a")







if __name__ == '__main__':
    main()


