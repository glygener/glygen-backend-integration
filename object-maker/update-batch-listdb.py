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
import csvutil
import libgly
import subprocess

def get_site_filter_code(obj):

    code_list = []

    #grp_idx=0 by_amino_acid
    c_list = []
    site_seq = obj["site_seq"]
    aa_three = obj["residue"] if "residue" in obj else ""
    aa_one =  aa_dict["one"][aa_three] if aa_three in aa_dict["one"] else "x"
    tmp_list = []
    if aa_one != "":
        tmp_list.append(aa_one)
    if len(site_seq) == 1 and site_seq not in tmp_list:
        tmp_list.append(site_seq)

    val_list = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    for val in val_list:
        c_list.append("1" if val in tmp_list else "0")
    code_list.append("".join(c_list))
    
    #grp_idx=1 by_glycosylation_type
    c_list = []
    tmp_list = obj["glycosylation_type"].replace(" ", "").split(";")
    val_list = ["n-linked|complex","n-linked|high mannose","n-linked|hybrid","n-linked|other"]
    val_list += ["o-linked|o-fucosylation", "o-linked|o-galactosylation","o-linked|o-galnacylation"]
    val_list += ["o-linked|o-glcnacylation","o-linked|o-glucosylation","o-linked|o-mannosylation"]
    val_list += ["o-linked|other","c-linked|c-mannosylation","c-linked|other","s-linked|other"]
    for val in sorted(val_list):
        c_list.append("1" if val in tmp_list else "0")
    code_list.append("".join(c_list))

    #grp_idx=2 by_organism filter
    c_list = []
    tmp_list = obj["organism"].replace(" ", "").split(";")
    for s in species_list:
        c_list.append("1" if s in tmp_list else "0")
    code_list.append("".join(c_list))


    #grp_idx=3 by_site_type
    c_list = []
    for f in ["glycation", "glycosylation", "mutagenesis", "phosphorylation", "snv"]:
        c_list.append("1" if obj[f] == "yes" else "0")
    code_list.append("".join(c_list))


    #grp_idx=4 by_snv_type
    c_list = []
    tmp_list = obj["snv_type"].replace(" ", "").split(";")
    for val in ["germline", "somatic"]:
        c_list.append("1" if val in tmp_list else "0")
    code_list.append("".join(c_list))


    #grp_idx=5 by_evidence_type filters
    c_list = []
    reported_flag = "1" if obj["reported_glycosite_flag"] == True else "0"
    mined_flag = "1" if obj["mined_glycosite_flag"] == True else "0"
    predicted_flag = "1" if obj["predicted_glycosite_flag"] == True else "0"
    c_list = [reported_flag, mined_flag, predicted_flag]
    code_list.append("".join(c_list))



    
    return ".".join(code_list)





def get_disease_filter_code(obj):
    
    code_list = []
    #grp_idx=0 by_biomarker_type
    c_list = []
    f_list = ["diagnostic", "monitoring", "predictive", "prognostic", "risk", "response", "safety"]
    type_list = obj["biomarker_type"].replace(" ", "").split(";")
    for f in f_list:
        c_list.append("1" if f in type_list else "0")
    code_list.append("".join(c_list))

    #grp_idx=1 by_data
    c_list = []
    f_list = ["protein_count", "glycan_count", "biomarker_count"]
    for f in f_list:
        c_list.append("1" if obj[f] > 0 else "0")
    code_list.append("".join(c_list))


    #grp_idx=3 by_organism filter
    c_list = []
    #tmp_list = obj["species_list"].replace(" ", "").split(";")
    tmp_list = obj["species_list"].split(";")
    for s in species_list:
        f = "1" if s in tmp_list else "0"
        c_list.append(f)
    code_list.append("".join(c_list))


    return ".".join(code_list)




def get_glycan_filter_code(obj):
    code_list = []

    #grp_idx=0 by_biomarker_type
    c_list = []
    f_list = ["diagnostic", "monitoring", "predictive", "prognostic", "risk", "response", "safety"]
    type_list = obj["biomarker_type"].replace(" ", "").split(";")
    for f in f_list:
        c_list.append("1" if f in type_list else "0")
    code_list.append("".join(c_list))

    #grp_idx=1 by_data
    c_list = []
    f_list = ["number_enzymes", "glycoprotein_count", "interactions_count", "publication_count"]
    for f in f_list:
        c_list.append("1" if obj[f] > 0 else "0")
    code_list.append("".join(c_list))


    #grp_idx=2 by_glycan_type
    c_list = []
    seen_val = {}
    for val in obj["classification"].replace(" ", "").split(";"):
        val = val.split("/")[0]
        seen_val[val] = True
    for val in ["N-linked","O-linked","Glycosphingolipid", "Other"]:
        c_list.append("1" if val in seen_val else "0")
    code_list.append("".join(c_list))

    #grp_idx=3 by_mass filters
    range_list = [
        {"min":0.0, "max":1000.0},
        {"min":1000.0, "max":2000.0},
        {"min":2000.0, "max":3000.0},
        {"min":3000.0, "max":4000.0},
        {"min":4000.0, "max":5000.0},
        {"min":5000.0, "max":100000000.0}
    ]
    c_list = []
    mass = obj["mass"]  if "mass" in obj else -1.0
    for o in range_list:
        c_list.append("1" if mass > o["min"] and mass < o["max"] else "0")
    c_list.append("1" if mass == -1.0 else "0")
    code_list.append("".join(c_list)) 

    #grp_idx=4 by_monosaccharide
    c_list = []
    tmp_list = obj["composition"].split(";")
    seen_val = {}
    for val in tmp_list:
        val = val.strip().split(" ")[0]
        seen_val[val] = True
    for val in ["hex","hexnac","dhex","hexa","neuac","neugc","hexn","p","s","other"]:
        c_list.append("1" if val in seen_val else "0")
    code_list.append("".join(c_list))    


    #grp_idx=5 by_organism filter
    c_list = []
    #tmp_list = obj["species_list"].replace(" ", "").split(";")
    tmp_list = obj["species_list"].split(";")
    
    for s in species_list:
        f = "1" if s in tmp_list else "0"
        c_list.append(f)
    code_list.append("".join(c_list))

 
    #grp_idx=6 by_sequence_details
    c_list = []
    tmp_list = obj["sequence_details"].strip().split(";")
    for val in ["Fully defined", "Topology", "Composition", "BaseComposition", "Incomplete"]:
        c_list.append("1" if val in tmp_list else "0")
    code_list.append("".join(c_list))
    


    return ".".join(code_list)



def get_biomarker_filter_code(obj):

    code_list = []

    #grp_idx=0 by_assessed_entity_type
    c_list = []
    f_list = ["glycan", "protein"]
    type_list = obj["assessed_entity_type"].replace(" ", "").split(";")
    for f in f_list:
        c_list.append("1" if f.lower() in type_list else "0")
    code_list.append("".join(c_list))
   
 
    #grp_idx=1 by_biomarker_type
    c_list = []
    f_list = ["diagnostic", "monitoring", "predictive", "prognostic", "risk", "response", "safety"]
    type_list = obj["best_biomarker_role"].replace(" ", "").split(";")
    for f in f_list:
        c_list.append("1" if f in type_list else "0")
    code_list.append("".join(c_list))

    return ".".join(code_list)


def get_protein_filter_code(obj):
 
    code_list = [] 
   
    #grp_idx=0 by_biomarker_type
    c_list = []
    f_list = ["diagnostic", "monitoring", "predictive", "prognostic", "risk", "response", "safety"]

    type_list = obj["biomarker_type"].replace(" ", "").split(";")
    for f in f_list:
        c_list.append("1" if f in type_list else "0")
    code_list.append("".join(c_list))



    #grp_idx=1 by_data
    c_list = []
    f_list = ["reported_snv","reported_mutagensis","reported_glycation","reported_interactions","publication_count"]
    f_list += ["disease_count", "pathway_count"]
    for f in f_list:
        c_list.append("1" if obj[f] > 0 else "0")
    code_list.append("".join(c_list))

    #grp_idx=2 by_mass filters
    range_list = [
        {"min":0.0, "max":100000.0},
        {"min":100000.0, "max":300000.0},
        {"min":300000.0, "max":500000.0},
        {"min":500000.0, "max":750000.0},
        {"min":750000.0, "max":100000000.0}
    ]
    c_list = []
    for o in range_list:
        c_list.append("1" if obj["mass"] > o["min"] and obj["mass"] < o["max"] else "0")
    c_list.append("1" if obj["mass"] == -1.0 else "0")
    code_list.append("".join(c_list))


    #grp_idx=3 by_organism filter
    c_list = []
    for s in species_list:
        c_list.append("0" if s != obj["organism"] else "1")
    code_list.append("".join(c_list))


    #grp_idx=4 by_ptm filters
    c_list = []
    f_list = ["total_reported_n_glycosites","total_reported_o_glycosites", "predicted_n_glycosites","predicted_o_glycosites"]
    f_list += ["reported_phosphosites"]
    for f in f_list:
        c_list.append("1" if obj[f] > 0 else "0")
    code_list.append("".join(c_list))



    #grp_idx=5 by_evidence_type filters
    c_list = []
    reported_flag = "1" if obj["reported_n_glycosites"] + obj["reported_o_glycosites"] > 0 else "0"
    mined_flag = "1" if obj["mined_glycosites"] > 0 else "0"
    predicted_flag = "1" if obj["predicted_glycosites"] > 0 else "0"
    c_list = [reported_flag, mined_flag, predicted_flag]
    code_list.append("".join(c_list))




    return ".".join(code_list)






def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " )
    parser.add_option("-s","--start",action="store",dest="start",help="")
    parser.add_option("-e","--end",action="store",dest="end",help="")

    (options,args) = parser.parse_args()
    for file in ([options.start, options.end]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    start = int(options.start)
    end = int(options.end)

    global config_obj
    global path_obj
    global species_list
    global data_dir
    global main_dict
    global aa_dict
    global url_map

    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"

    #DEBUG = False
    DEBUG = True

    url_map = {}
    in_file = "generated/misc/list_init_urls.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        record_type = row[f_list.index("record_type")]
        field_id = row[f_list.index("field_id")]
        url_template = row[f_list.index("url_template")]
        if record_type not in url_map:
            url_map[record_type] = {}
        url_map[record_type][field_id] = url_template



    pattern_dict = {"protein":"*", "glycan":"*", "site":"*", "biomarker":"*", "motif":"*", "disease":"*"}
    if DEBUG:
        #pattern_dict = {"protein":"P41222*"}
        #pattern_dict = {"site":"P02766*"}
        #pattern_dict = {"site":"P41222*"}
        #pattern_dict = {"glycan":"*"}
        pattern_dict = {"disease":"*"}


    file_obj_list = []
    for record_type in sorted(list(pattern_dict.keys())):
        #glob_str = "jsondb/jumbodb/%sdb/%s.json" % (record_type, pattern_dict[record_type])
        glob_str = "jsondb/%sdb/%s.json" % (record_type, pattern_dict[record_type])
        file_list = glob.glob(glob_str)
        for json_file in sorted(file_list):
            file_obj_list.append({"file":json_file, "record_type":record_type})

    #print (len(file_obj_list))
    #exit()
    
    log_file = "logs/update-listdb.%s-%s.log" % (start, end)
    msg = "update-listdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")

    aa_dict = csvutil.load_aa_format_dict()
    species_list = []
    species_obj = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)
    for tax_id in species_obj:
        if species_obj[tax_id]["is_reference"] == "yes":
            name = species_obj[tax_id]["glygen_name"]
            if name not in species_list:
                species_list.append(name)
    species_list = sorted(species_list)
    total = len(file_obj_list[start-1:end])
    end = len(file_obj_list) if end > len(file_obj_list) else end


    record_count_dict = {}
    for obj in file_obj_list[start-1:end]:
        record_type = obj["record_type"]
        in_file = obj["file"]
        if record_type not in record_count_dict:
            record_count_dict[record_type] = 0
        r = record_count_dict[record_type]
        list_file = in_file.replace(record_type, "list").replace("jumbodb/", "")
        if os.path.isfile(list_file) == False:
            msg = "update-listdb: list_file not found for %s" % (in_file)
            csvutil.write_log_msg(log_file, msg, "a")
            continue
        doc = json.loads(open(list_file,"r").read())
        if r > 0 and r%1000 == 0:
            msg = "update-listdb: created %s/%s %s listdb records" % (r,total, record_type)
            csvutil.write_log_msg(log_file, msg, "a")
        doc = json.loads(open(list_file,"r").read())
        if record_type == "protein":
            doc["filter_code"] = get_protein_filter_code(doc)
        elif record_type in "glycan":
            doc["filter_code"] = get_glycan_filter_code(doc)
        elif record_type == "motif":
            doc["filter_code"] = "000.000.000"
        elif record_type == "site":
            doc["filter_code"] = get_site_filter_code(doc)
        elif record_type == "biomarker":
            doc["filter_code"] = get_biomarker_filter_code(doc)
        elif record_type == "disease":
            doc["filter_code"] = get_disease_filter_code(doc)
        if record_type in url_map:
            for field_id in url_map[record_type]:
                if field_id in doc:
                    new_field_id = field_id + "_url"
                    if new_field_id in doc:
                        doc.pop(new_field_id)
                    field_value = doc[field_id]
                    if type(field_value) is str:
                        oo_list = []
                        if field_value.strip() != "":
                            for val in field_value.strip().split(";"):
                                url = url_map[record_type][field_id] % (val)
                                oo_list.append({"id":val, "url":url})
                        #print (field_id, oo_list)
                        doc[field_id] = oo_list
        with open(list_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        record_count_dict[record_type] += 1  
   

    for record_type in record_count_dict:
        r = record_count_dict[record_type]
        msg = "update-listdb: final created: %s/%s %s listdb records" % (r, total, record_type)
        csvutil.write_log_msg(log_file, msg, "a")



if __name__ == '__main__':
    main()

