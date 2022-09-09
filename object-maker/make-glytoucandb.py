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
    
    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]


    in_file = "reviewed/glycan_glytoucanidlist.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    out_obj = {"aclist":[]}
    for row in data_frame["data"]:
        out_obj["aclist"].append(row[0])

    out_file = path_obj["jsondbpath"] + "/glytoucandb/glytoucandb.json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(out_obj, indent=4)))
    print ("make-glytoucandb: final created: %s glytoucandb objects" % (1))


if __name__ == '__main__':
    main()

