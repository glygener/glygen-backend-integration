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

    
    species_obj = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)



    record_count = 0
    seen_species = {}
    file_list = glob.glob("jsondb/proteindb/*.json")
    file_list += glob.glob("jsondb/glycandb/*.json")
    for in_file in file_list:
        doc = json.loads(open(in_file,"r").read())
        for o in doc["species"]:
            combo_id = o["taxid"]
            tax_id = str(combo_id)
            if tax_id in species_obj:
                if species_obj[tax_id]["is_reference"] == "no":
                    continue
            if combo_id not in seen_species:
                seen_species[combo_id] = True
                out_file = "jsondb/speciesdb/%s.json" % (combo_id)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(o, indent=4)))
                    record_count += 1 
    
    log_file = "logs/make-speciesdb.log"
    msg = "make-speciesdb: final created: %s species objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "w")





if __name__ == '__main__':
    main()

