#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import requests
import subprocess

def main():

    global data_path
    global data_root
    global data_rel

    global wrk_dir
    wrk_dir = "/home/rykahsay/glygen-backend-integration/object-maker"
    file_list = glob.glob("%s/html/*.html" % (wrk_dir))
    for in_file in file_list:
        page_id = in_file.split("/")[-1].split(".")[0]
        cn = open(in_file, "r").read();
        doc = {"pageid":page_id, "cn":cn}
        out_file = wrk_dir + "/jsondb/htmldb/%s.json" % (page_id)
        with open(out_file,"w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        print ("created jsondb/htmldb/%s.json" % (page_id))
    
    cmd = "chmod -R 775 " + wrk_dir + "/jsondb/htmldb"
    x, y = subprocess.getstatusoutput(cmd)

if __name__ == '__main__':
    main()


