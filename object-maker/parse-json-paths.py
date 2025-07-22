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



def parse_doc(in_obj, seen, path):
    

    if type(in_obj) is dict:
        for k in in_obj:
            p = path + "." + k if path != "" else k
            parse_doc(in_obj[k], seen, p)
    elif type(in_obj) is list:
        for k in range(0, len(in_obj)):
            p = path + "." + str(k) if path != "" else str(k)
            parse_doc(in_obj[k], seen, p)
    elif type(in_obj) in [int, float, str]:
        seen[path] = str(in_obj)
    

    return





def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog ")
    parser.add_option("-r","--recordtype",action="store",dest="recordtype")

    (options,args) = parser.parse_args()
    for file in ([options.recordtype]):
        if not (file):
            parser.print_help()
            sys.exit(0)



    record_type = options.recordtype
    log_file = "logs/value_stat_%s.log" % (record_type)
    out_file = "logs/final_value_stat_%s" % (record_type)

    file_list = glob.glob("jsondb/%sdb/*.json" % (record_type))
    stat_dict = {}
    FL = open(log_file, "w")
    n = 0
    for in_file in file_list:

        seen = {}
        doc = json.loads(open(in_file, "r").read())
        parse_doc(doc, seen, "")
        tmp_dict = {}
        for path in seen:
            tmp_list = []
            for k in path.split("."):
                tmp_list.append("*" if k.isdigit() else k)
            p = ".".join(tmp_list)
            tmp_dict[p] = seen[path]
        for p in tmp_dict:
            if p not in stat_dict:
                stat_dict[p] = 0
            stat_dict[p] += 1
        n += 1
        if n%1000 == 0:
            FL.write("... processed %s %s objects" % (n, record_type))
    FL.close()
    
    with open(out_file, "w") as FW:
        for p in stat_dict:
            FW.write("%s / %s : %s\n" % (stat_dict[p], n, p))




if __name__ == '__main__':
    main()

