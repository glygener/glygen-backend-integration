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




def get_by_amino_acid(aa_three, site_seq):


    aa_one =  aa_dict["one"][aa_three] if aa_three in aa_dict["one"] else "x"
    tmp_list = []
    if aa_one != "":
        tmp_list.append(aa_one)
    if len(site_seq) == 1 and site_seq not in tmp_list:
        tmp_list.append(site_seq)
    val_list = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    c_list = []
    for val in val_list:
        c_list.append("1" if val in tmp_list else "0")
    return "".join(c_list)

def get_by_glycosylation_type(glycosylation_type):

    c_list = []
    tmp_list = glycosylation_type.replace(" ", "").split(";")
    val_list = ["n-linked|complex","n-linked|high mannose","n-linked|hybrid","n-linked|other"]
    val_list += ["o-linked|o-fucosylation", "o-linked|o-galactosylation","o-linked|o-galnacylation"]
    val_list += ["o-linked|o-glcnacylation","o-linked|o-glucosylation","o-linked|o-mannosylation"]
    val_list += ["o-linked|other","c-linked|c-mannosylation","c-linked|other","s-linked|other"]
    for val in val_list:
        c_list.append("1" if val in tmp_list else "0")
    return "".join(c_list)


def get_by_prd_tool(prd_tool):
    c_list = []
    val_list = ["isoglyp"]
    for val in val_list:
        c_list.append("1" if val == prd_tool.lower() else "0")
    return "".join(c_list)



def get_by_mining_tool(mining_tool):
    m_list = mining_tool.lower().split(";")
    c_list = []
    val_list = ["RLIMS-G", "GlycoSiteMiner"]
    for val in val_list:
        c_list.append("1" if val.lower() in m_list else "0")
    return "".join(c_list)



def get_code_list(canon, obj, code_dict):

    
    #glycan order:    by_biomarker_type by_data by_glycan_type by_mass by_monosaccharide by_organism by_sequence_details
    #protein order:   by_biomarker_type by_data by_mass by_organism by_ptm
    #site order:      by_amino_acid by_glycosylation_type by_organism by_site_type by_snv_type
    #site_w_glycan:   by_amino_acid by_biomarker_type by_glycosylation_type by_mass by_monosaccharide 
    #                 by_organism by_sequence_details by_site_type
    #site_reported:   by_amino_acid by_glycosylation_type by_site_type 
    #site_predicted:  by_amino_acid by_glycosylation_type by_prd_tool by_site_type

    tmpl_dict = {}
    tmpl_dict["glycan"] = "0000000.0000.0000.0000000.0000000000.00000000000000.00000".replace("0","x").split(".")
    tmpl_dict["site"] = "00000000000000000000.00000000000000.00000000000000.00000.00".replace("0","x").split(".")

    site_category_dict = obj["site_category_dict"]
    start_pos = obj["start_pos"] if "start_pos" in obj else "0"
    end_pos = obj["start_pos"] if "end_pos" in obj else "0"
    aa_three = obj["residue"] if "residue" in obj else ""
    site_seq = obj["site_seq"] if "site_seq" in obj else ""
    prd_tool = obj["prediction_tool"] if "prediction_tool" in obj else ""
    mining_tool = obj["mining_tool"] if "mining_tool" in obj else ""
    site_id = "%s.%s.%s" % (canon, start_pos, end_pos)
    obj["subtype"] = "other" if obj["subtype"] == "" else obj["subtype"]
    gly_type = "%s|%s" % (obj["type"].lower(), obj["subtype"].lower())
    site_code_list =  code_dict[site_id].split(".") if site_id in code_dict else tmpl_dict["site"]
    if "reported_with_glycan" in site_category_dict:
        glycan_id = obj["glytoucan_ac"]
        glycan_code_list =  code_dict[glycan_id].split(".") if glycan_id in code_dict else tmpl_dict["glycan"]
        code_list = ["x", glycan_code_list[0], "x",glycan_code_list[3], 
                    glycan_code_list[4],glycan_code_list[5], glycan_code_list[6], site_code_list[3]]
        aa_three = obj["residue"] if "residue" in obj else ""
        code_list[0] = get_by_amino_acid(aa_three, site_seq)
        code_list[2] = get_by_glycosylation_type(gly_type)
        obj["filter_code"] = ".".join(code_list)
        #print (site_id,glycan_id,obj["filter_code"]) 
    elif "reported" in site_category_dict or "automatic_literature_mining" in site_category_dict:
        code_list = ["x", "x", "x", "x"]
        code_list[0] = get_by_amino_acid(aa_three, site_seq)
        code_list[1] = get_by_glycosylation_type(gly_type)
        code_list[2] = get_by_mining_tool(mining_tool)
        code_list[3] = site_code_list[3]
        obj["filter_code"] = ".".join(code_list)
        print (site_id, mining_tool, obj["filter_code"]) 
    elif "predicted" in site_category_dict:
        code_list = ["x", "x", "x", "x"]
        code_list[0] = get_by_amino_acid(aa_three, site_seq)
        code_list[1] = get_by_glycosylation_type(gly_type)
        code_list[2] = get_by_prd_tool(prd_tool)
        code_list[3] = site_code_list[3]
        obj["filter_code"] = ".".join(code_list)
        #print (site_id, site_category, obj["filter_code"])

    return code_list






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


    global species_list
    global code_dict
    global aa_dict
    

    DEBUG = False
    #DEBUG = True

    pattern_dict = {"protein":"*"}
    if DEBUG:
        #pattern_dict = {"protein":"P20827*", "glycan":"G17689*", "site":"P14210*"}
        pattern_dict = {"protein":"P41222*"}


    
    log_file = "logs/update-proteindb.%s-%s.log" % (start, end)
    msg = "update-proteindb: started logging"
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


    code_dict = {}
    record_type_list = ["protein", "glycan", "site"]
    for record_type in record_type_list:
        file_list = glob.glob("jsondb/%sdb/*" % (record_type))
        idx = 0
        #Loading filter_codes from $record_type list objects
        for in_file in file_list:
            list_file = in_file.replace(record_type, "list")  
            if os.path.isfile(list_file) == False:
                msg = "update-proteindb: list_file not found for %s" % (in_file.split("/")[-1])
                csvutil.write_log_msg(log_file, msg, "a")
                continue
            if idx%10000 == 0:
                msg = "update-proteindb: loaded %s/%s %s filter codes" % (idx, len(file_list), record_type)
                csvutil.write_log_msg(log_file, msg, "a")
            doc = json.loads(open(list_file,"r").read())
            record_id = doc["record_id"]
            if "filter_code" in doc:
                code_dict[record_id] = doc["filter_code"]
            idx += 1
 
    record_type = "protein"
    file_list = glob.glob("jsondb/%sdb/%s" % (record_type, pattern_dict[record_type]))
    end = end if end < len(file_list) else len(file_list)


    total = len(file_list)
    record_count = 0
    for in_file in file_list[start-1:end]:
        doc = json.loads(open(in_file,"r").read())
        canon = doc["uniprot_canonical_ac"]
        for obj in doc["glycosylation"]:
            code_list = get_code_list(canon, obj, code_dict)
            obj["filter_code"] = ".".join(code_list)
        with open(in_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        record_count += 1
        batch_file_list = glob.glob("jsondb/batchdb/protein.%s.*.json" % (canon))
        for batch_file in batch_file_list:
            batch_doc = json.loads(open(batch_file,"r").read())
            if "glycosylation" in batch_doc["sections"]:
                for obj in batch_doc["sections"]["glycosylation"]:
                    code_list = get_code_list(canon, obj, code_dict)
                    obj["filter_code"] = ".".join(code_list)
                with open(batch_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(batch_doc, indent=4)))
                record_count += 1
        if record_count%1000 == 0:
            msg = "update-proteindb: updated: %s/%s records" % (record_count, total)
            csvutil.write_log_msg(log_file, msg, "a")
 
    msg = "update-proteindb: final updated: %s/%s records" % (record_count, total)
    csvutil.write_log_msg(log_file, msg, "a")


if __name__ == '__main__':
    main()

