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


def main():

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global misc_dir
    global main_dict


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"
    misc_dir = "generated/misc/"


    record_count = 0
    seen_species = {}
    file_list = glob.glob("jsondb/proteindb/*.json")
    file_list += glob.glob("jsondb/glycandb/*.json")
    for in_file in file_list:
        doc = json.loads(open(in_file,"r").read())
        for o in doc["species"]:
            combo_id = o["taxid"]
            if combo_id not in seen_species:
                seen_species[combo_id] = True
                out_file = path_obj["jsondbpath"] + "/speciesdb/%s.json" % (combo_id)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(o, indent=4)))
                    record_count += 1 
    
    print ("make-speciesdb: final created: %s species objects" % (record_count))




if __name__ == '__main__':
    main()

