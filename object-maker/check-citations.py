#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import requests
import subprocess
from Bio import SeqIO

import csvutil
import libgly
import biomarker_util




def main():
   
    seen_one = {}
    file_list = glob.glob("jsondb/publicationdb/*.json")
    for in_file in file_list:
        file_name = in_file.split("/")[-1].replace(".json", "")
        seen_one[file_name] = True

    seen_two = {}
    #file_list = glob.glob("reviewed/*.csv")
    tmp_dict = {}
    for line in open("tmp/JUNK", "r").read().split("\n"):
        file_name = line.split("|")[-1]
        if file_name != "":
            f = "reviewed/" + file_name
            tmp_dict[f] = True

    file_list = list(tmp_dict.keys())

    n = len(file_list)
    idx = 0
    for in_file in file_list:
        idx += 1
        if in_file.find(".stat.csv") != -1:
            continue
        file_name = in_file.split("/")[-1]
        if file_name in ["glycan_synthesized.csv"]:
            continue
        data_frame = {}
        csvutil.load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        if "xref_key" not in f_list:
            continue
        for row in data_frame["data"]:
            xref_key, xref_id = row[f_list.index("xref_key")], row[f_list.index("xref_id")]
            k = xref_key.split("_")[-1]
            if k not in ["pubmed", "doi"]:
                continue
            if k == "doi":
                xref_id = xref_id.replace("/", "_").lower()
            c = "%s.%s" % (k, xref_id)
            cc = "%s|%s|%s" % (k,xref_id,file_name)
            if c not in seen_one and cc not in seen_two:
                seen_two[cc] = True
                print (cc)
    exit()

    for cc in seen_two:
        db, xref_id, file_name = cc.split("|")
        cmd = "grep \"%s\" reviewed/%s |head -1  " % (xref_id, file_name)
        x = subprocess.getoutput(cmd).strip().split(",")[0]
        print (x, cc)


if __name__ == '__main__':
    main()



