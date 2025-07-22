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


def sort_release_list(tmp_list, reversed_flag):

    factor_list = [100000000, 1000, 1]
    rel_dict = {}
    for rel in tmp_list:
        parts = rel.split(".") if rel.find(".") != -1 else rel.split("_")
        ordr = 0
        for i in range(0,len(parts)):
            p = int(parts[i]) if parts[i] != "x" else 100
            ordr += factor_list[i]*p
        rel_dict[ordr] = rel
    
    release_list = []

    for ordr in sorted(rel_dict, reverse=reversed_flag):
        release_list.append(rel_dict[ordr])

    return release_list






def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " )
    parser.add_option("-v","--ver",action="store",dest="ver",help="2.0.2")

    (options,args) = parser.parse_args()
    for file in ([options.ver]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    global wrk_dir
    
    ver = options.ver

    dir_list = glob.glob("jsondb/*")
    for d_new in dir_list:
        db = d_new.split("/")[-1]
        d_old = "/data/shared/glygen/releases/data/v-%s/jsondb/%s/" % (ver, db)
        
        n_new = len(glob.glob(d_new + "/*"))
        n_old = 0
        if os.path.isdir(d_old):
            n_old = len(glob.glob(d_old + "/*"))
        print ("%s: %s --> %s [%s]" % (db, n_old, n_new, n_new - n_old))



if __name__ == '__main__':
    main()

