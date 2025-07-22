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
import csvutil




def get_enum_dict():

    doc_dict = {}
    for record_type in ["species", "glycan", "motif", "protein"]:
        doc_dict[record_type] = []
        file_list = glob.glob("jsondb/%sdb/*.json" % (record_type))
        for in_file in file_list:
            doc = json.loads(open(in_file, "r").read())
            doc_dict[record_type].append(doc)

    enum_dict = {
        "glycan":{},
        "motif":{},
        "protein":{},
        "enzyme":{},
        "site":{},
        "species":{},
        "gene":{},
        "disease":{}
    }
 
    # glycan_classification
    record_type = "glycan"
    coll = "c_glycan"
    f_one = "classification.type.name"
    f_two = "classification.subtype.name"
    f_three = "crossref.database"
    enum_dict[record_type][f_one] = []
    enum_dict[record_type][f_two] = []
    enum_dict[record_type][f_three] = []
    for doc in doc_dict["glycan"]:
        for o in doc["classification"]:
            if o["type"]["name"] not in enum_dict[record_type][f_one]:
                enum_dict[record_type][f_one].append(o["type"]["name"])
            if o["subtype"]["name"] not in enum_dict[record_type][f_two]:
                enum_dict[record_type][f_two].append(o["subtype"]["name"])
        for o in doc["crossref"]:
            val = o["database"]
            if val not in enum_dict[record_type][f_three]:
                enum_dict[record_type][f_three].append(val)


    enum_dict["protein"][f_three] = []
    enum_dict["enzyme"][f_three] = []
     
    for doc in doc_dict["protein"]:
        for o in doc["crossref"]:
            val = o["database"]
            record_type = "protein"
            if val not in enum_dict[record_type][f_three]:
                enum_dict[record_type][f_three].append(val)
            if "enzyme" in doc["keywords"]:
                record_type = "enzyme"
                if val not in enum_dict[record_type][f_three]:
                    enum_dict[record_type][f_three].append(val)

    # motif aglycon, reducing_end, alignment_method
    record_type = "motif"
    coll = "c_motif"
    f_one, f_two, f_three = "aglycon", "reducing_end", "alignment_method"
    enum_dict[record_type][f_one] = []
    enum_dict[record_type][f_two] = []
    enum_dict[record_type][f_three] = []
    for doc in doc_dict["motif"]:
        for f in [f_one, f_two, f_three]:
            val = doc[f].strip() if f in doc else ""
            if val != "" and val not in enum_dict[record_type][f]:
                enum_dict[record_type][f].append(val)

    site_type_list = [
        "glycosylation_flag", "snv_flag", "phosphorylation_flag", "glycation_flag",
        "mutagenesis_flag",
        "glycosylation", "snv", "phosphorylation", "glycation",
        "mutagenesis"
    ]

    for f in site_type_list:
        enum_dict["site"][f] = ["true", "false"]
    for f in ["fully_determined"]:
        enum_dict["glycan"][f] = ["yes", "no"]
    for f in ["neighbors.direction"]:
        enum_dict["site"][f] = ["upstream", "downstream"]
    for f in ["neighbors.categories"]:
        enum_dict["site"][f] = site_type_list

    enum_dict["species"]["taxid"], enum_dict["species"]["name"] = [], []
    for doc in doc_dict["species"]:
        if "glygen_name" in doc:
            enum_dict["species"]["taxid"].append(doc["taxid"])
            enum_dict["species"]["name"].append(doc["glygen_name"])


    return enum_dict





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
            path_info[record_type][path] = {"label":label, "description":desc, "order":order}


    return


def main():


    path_info = {}
    load_path_info(path_info)
    
    enum_dict = get_enum_dict()

    record_count = 0
    obj_dict = json.loads(open("generated/misc/prop_master.json", "r").read())
    for record_type in obj_dict:
        for doc in obj_dict[record_type]:
            label, desc, order = "", "", 0
            path = doc["path"]
            if path in path_info[record_type]:
                if "label" in path_info[record_type][path]:
                    label = path_info[record_type][path]["label"].strip()
                if "description" in path_info[record_type][path]:
                    desc = path_info[record_type][path]["description"].strip()
                if "order" in path_info[record_type][path]:
                    if str(path_info[record_type][path]["order"]).strip().isdigit() == True:
                        order = int(str(path_info[record_type][path]["order"]).strip()) 
            doc["label"] = label if label != "" else doc["label"]
            doc["description"] = desc if desc != "" else doc["description"]
            doc["enum"] = []
            if record_type in enum_dict:
                if path in enum_dict[record_type]:
                    doc["enum"] = enum_dict[record_type][path]
            doc["order"] = order if order != 0 else doc["order"]
            json_file = "jsondb/pathdb/%s.%s.json" % (record_type, path)
            with open(json_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            record_count += 1

    log_file = "logs/make-pathdb.log"
    msg = "make-pathdb: ... final created: %s path objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "w")



if __name__ == '__main__':
    main()

