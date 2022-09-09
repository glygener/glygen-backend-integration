#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import subprocess
from optparse import OptionParser

__version__="1.0"
__status__ = "Dev"






def main():


    global wrk_dir

    wrk_dir = "/home/rykahsay/glygen-backend-integration/object-maker"
    file_list = glob.glob(wrk_dir + "/jsondb/bcodb/*.json")
    for in_file in sorted(file_list):
        doc = json.loads(open(in_file, "r").read())
        bco_id = doc["object_id"].split("/")[-2]
        for obj in doc["io_domain"]["output_subdomain"]:
            file_name = obj["uri"]["uri"].split("/")[-1]
            if file_name.find(".stat.csv") == -1 and file_name.split(".")[-1] != ".log":
                print (bco_id, file_name)






if __name__ == '__main__':
    main()


