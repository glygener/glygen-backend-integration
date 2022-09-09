#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import requests
import subprocess
import pymongo
from optparse import OptionParser
import libgly
from Bio import SeqIO


__version__="1.0"
__status__ = "Dev"





def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-v","--ver",action="store",dest="ver",help="New version")

    (options,args) = parser.parse_args()
    for key in ([options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    global wrk_dir 
    old_rel = options.ver
    wrk_dir = "/home/rykahsay/glygen-backend-integration/object-maker"
    old_path_list = glob.glob("/data/shared/glygen/releases/data/v-%s/reviewed/*.*" % (old_rel))
    new_path_list = glob.glob(wrk_dir + "/reviewed/*.*")

    old_file_list = []
    for path in old_path_list:
        ignore =  False
        for k in [".stat.csv", "GLY_", "release-notes.txt"]:
            if path.find(k) != -1:
                ignore = True
        if ignore == True:
            continue
        file_name = path.split("/")[-1]
        old_file_list.append(file_name)

    new_file_list = []
    for path in new_path_list:
        ignore =  False
        for k in [".stat.csv", "GLY_", "release-notes.txt"]:
            if path.find(k) != -1:
                ignore = True
        if ignore == True:
            continue
        file_name = path.split("/")[-1]
        new_file_list.append(file_name)

    for file_name in old_file_list:
        if file_name not in new_file_list:
            print ("missing file:", file_name)

    for file_name in new_file_list:
        if file_name not in old_file_list:
            print ("new file:", file_name)



if __name__ == '__main__':
    main()



