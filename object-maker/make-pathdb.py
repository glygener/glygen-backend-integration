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
import collections


import libgly


def load_path_info(path_info):

    file_list = glob.glob("generated/misc/*_propertyinfo.csv")
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

    return


def main():


    path_info = {}
    load_path_info(path_info)

    record_count = 0
    obj_dict = json.loads(open("generated/misc/prop_master.json", "r").read())
    for record_type in obj_dict:
        for doc in obj_dict[record_type]:
            label, desc, order = "", "", 0
            path = doc["path"]
            if path in path_info:
                if "label" in path_info[path]:
                    label = path_info[path]["label"].strip()
                if "description" in path_info[path]:
                    desc = path_info[path]["description"].strip()
                if "order" in path_info[path]:
                    if str(path_info[path]["order"]).strip().isdigit() == True:
                        order = int(str(path_info[path]["order"]).strip()) 
            doc["label"] = label if label != "" else doc["label"]
            doc["description"] = desc if desc != "" else doc["description"]

            doc["order"] = order if order != 0 else doc["order"]
            json_file = "jsondb/pathdb/%s.%s.json" % (record_type, path)
            with open(json_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            record_count += 1

    print ("make-pathdb: ... final created: %s path objects" % (record_count))


if __name__ == '__main__':
    main()

