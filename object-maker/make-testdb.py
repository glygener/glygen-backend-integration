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

import datetime
import pytz

import libgly
import csvutil

import subprocess
import hashlib


def parse_doc(in_obj, path_dict, path):
    

    if type(in_obj) is dict:
        for k in in_obj:
            p = path + "." + k if path != "" else k
            parse_doc(in_obj[k], path_dict, p)
    elif type(in_obj) is list:
        for k in range(0, len(in_obj)):
            p = path + "." + str(k) if path != "" else str(k)
            parse_doc(in_obj[k], path_dict, p)
    elif type(in_obj) is str:
        path_dict[path] = "*"
    elif type(in_obj) is int:
        path_dict[path] = "0"
    elif type(in_obj) is float:
        path_dict[path] = "0.0"

    return





def main():

    batch_size = 1000
    record_type_list = ["protein", "glycan"]
    #record_type_list = ["protein"]

    log_file = "logs/make-testdb.log" 
    with open(log_file, "w") as FL:
        FL.write("Started logging\n")

    for record_type in record_type_list:
        file_list = glob.glob("jsondb/%sdb/*.json" % (record_type))
        stat_dict = {}
        n_records, n_pass = 0, 0
        seen = {}
        bucket_list = []
        for in_file in file_list:
            doc = json.loads(open(in_file, "r").read())
            main_id = doc["uniprot_canonical_ac"] if record_type == "protein" else doc["glytoucan_ac"]
            path_dict = {}
            parse_doc(doc, path_dict, "")
            tmp_dict = {}
            for path in path_dict:
                tmp_list = []
                for k in path.split("."):
                    tmp_list.append("*" if k.isdigit() else k)
                new_path = ".".join(tmp_list)
                tmp_dict[new_path] = path_dict[path]
            hash_str = json.dumps(tmp_dict)
            hash_obj = hashlib.md5(hash_str.encode('utf-8'))
            digest = hash_obj.hexdigest()
            if digest not in seen:
                bucket_list.append(main_id)
                seen[digest] = True
                n_pass += 1
            if n_records > 0 and n_records%1000 == 0:
                with open(log_file, "a") as FL:
                    FL.write("... processed %s %s objects, %s passed\n" % (n_records, record_type, n_pass))
            n_records += 1

        start_idx = 0
        batch_idx = 1
        n = len(bucket_list)
        while start_idx < n:
            end_idx = start_idx + batch_size
            end_idx = n if end_idx > n else end_idx
            tmp_list = bucket_list[start_idx:end_idx]
            out_file = "jsondb/testdb/%s.batch.%s.json" % (record_type, batch_idx)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(tmp_list, indent=4)))
            start_idx += batch_size
            batch_idx += 1





if __name__ == '__main__':
    main()

