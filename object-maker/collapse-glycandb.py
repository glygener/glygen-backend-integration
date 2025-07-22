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


import libgly
import csvutil
import section_stats


def collapse_objects(obj_list):

    seen = {}
    obj_dict = {}
    for obj in obj_list:
        canon = obj["uniprot_canonical_ac"] if "uniprot_canonical_ac" in obj else ""
        aa_pos = obj["start_pos"] if "start_pos" in obj else ""
        tissue_obj = obj["tissue"]
        cell_line_obj = obj["cell_line"]
        sample_src_id = cell_line_obj["id"] if tissue_obj == {} else tissue_obj["id"]
        combo_id_one = "%s|%s|%s" % (canon,aa_pos,sample_src_id)
        if combo_id_one not in seen:
            if "abundance" in obj:
                obj.pop("abundance")
            seen[combo_id_one] = {}
            obj_dict[combo_id_one] = obj
        for ev_obj in obj["evidence"]:
            combo_id_two = json.dumps(ev_obj)
            seen[combo_id_one][combo_id_two] = True
    
    new_obj_list = []
    for combo_id_one in seen:
        obj = obj_dict[combo_id_one]
        obj["evidence"] = []
        for combo_id_two in seen[combo_id_one]:
            obj["evidence"].append(json.loads(combo_id_two))
        new_obj_list.append(obj)
    
    return new_obj_list



def main():

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global main_dict


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"


    sec_list = config_obj["sitesections"]

    protein_obj_dict = {}
    record_count = 0

    DEBUG = False
    #DEBUG = True

    file_list = glob.glob("jsondb/glycandb/*.json")
    if DEBUG:
        file_list = glob.glob("jsondb/glycandb/G59655SA.json")
    for in_file in file_list:
        doc = json.loads(open(in_file,"r").read())
        if len(doc["expression"]) == 0:
            continue
        doc["expression"] = collapse_objects(doc["expression"])
        section_stats.get_sec_stats(doc, "glycan")
        with open(in_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))

    file_list = glob.glob("jsondb/batchdb/glycan.*.json")
    for in_file in file_list:
        doc = json.loads(open(in_file,"r").read())
        if "expression" in doc["sections"]:
            doc["sections"]["expression"] = collapse_objects(doc["sections"]["expression"])
            with open(in_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))



if __name__ == '__main__':
    main()

