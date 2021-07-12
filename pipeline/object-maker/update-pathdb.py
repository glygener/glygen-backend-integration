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
import collections


sys.path.append('../../glytools/')
import libgly



def load_properity_lineage(in_obj, in_key, seen):

    if type(in_obj) in [dict, collections.OrderedDict]:
        for k in in_obj:
            new_key = in_key + "." + k if in_key != "" else k
            load_properity_lineage(in_obj[k], new_key, seen)
    elif type(in_obj) is list:
        for o in in_obj:
            load_properity_lineage(o, in_key, seen)
    elif type(in_obj) in [unicode, int, float]:
        value_type = str(type(in_obj)).replace("<type ", "").replace("'", "").replace(">", "")
        in_key += " %s" % (value_type)
        if in_key not in seen:
            seen[in_key] = True

    return




def main():

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global misc_dir
    global main_dict


    config_file = "../../conf/config-1.1.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"
    misc_dir = "generated/misc/"

    file_list = glob.glob("generated/misc/*_propertyinfo.csv")
    path_info = {}
    for in_file in file_list:
        record_type = in_file.split("/")[-1].split("_")[0]
        line_list = open(in_file, "r").read().split("\n")
        for line in line_list[1:]:
            w_list = line.strip().split(",")
            if len(w_list) < 2:
                continue
            if record_type not in path_info:
                path_info[record_type] = {}
            path, label = w_list[0].split(" ")[0], w_list[1]
            desc = w_list[2] if len(w_list) > 2 else "xxx"
            order = w_list[3] if len(w_list) > 3 else 1000
            path_info[path] = {"label":label, "description":desc, "order":order}

    file_list = glob.glob("jsondb/pathdb/*.json") 
    for json_file in file_list:
        doc = json.loads(open(json_file,"r").read())
        record_type = json_file.split("/")[-1].split(".")[0]
         
        path = doc["path"].split(" ")[0]
        label, desc, order = "xxx", "xxx", 1000
        if path in path_info:
            if "label" in path_info[path]:
                label = path_info[path]["label"].strip()
            if "description" in path_info[path]:
                desc = path_info[path]["description"].strip()
            if "order" in path_info[path]:
                if str(path_info[path]["order"]).strip().isdigit() == True:
                    order = int(str(path_info[path]["order"]).strip()) 
        
        doc["path"] = path
        doc["label"] = label
        doc["description"] = desc
        doc["order"] = order
        with open(json_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))




if __name__ == '__main__':
    main()

