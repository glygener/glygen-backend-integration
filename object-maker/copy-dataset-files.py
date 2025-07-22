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



def get_master_file_list():

    file_name_list = []
    ds_obj_list = json.loads(open(wrk_dir + "/generated/misc/dataset-masterlist.json", "r").read())
    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if ds_name in ["homolog_alignments", "isoform_alignments"]:
            continue
        if obj["categories"]["species"] == []:
            file_name_list.append("%s_%s.%s" % (mol, ds_name, ds_format))
        else:
            sp_list_one = sorted(obj["categories"]["species"])
            for species in sp_list_one:
                if species not in obj["integration_status"]["excludelist"]:
                    file_name_list.append("%s_%s_%s.%s" % (species, mol, ds_name, ds_format))
    
    return file_name_list





def main():

    global wrk_dir
    global field_dict
    global io_dict 


    generated_dir = "/data/projects/glygen/generated/"
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    reviewed_dir = wrk_dir + "/reviewed/"
    unreviewed_dir = wrk_dir + "/unreviewed/"

    file_list = get_master_file_list()

    path_list = []
    missing_files = []
    for out_file_name in file_list:
        path = unreviewed_dir + out_file_name
        if  os.path.isfile(path) == False:
            missing_files.append(path)
        else:
            path_list.append(path)
   

    if missing_files != []:
        for path in missing_files:
            print (path, "is missing")
    else:
        cmd = "rm -f " + reviewed_dir + "/*"
        x, y = subprocess.getstatusoutput(cmd)
        for path in path_list:
            cmd = "cp " + path + " "  + reviewed_dir
            x, y = subprocess.getstatusoutput(cmd)

    cmd = "chmod -R 755 " + reviewed_dir
    x, y = subprocess.getstatusoutput(cmd)


if __name__ == '__main__':
    main()



