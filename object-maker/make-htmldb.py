#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import requests
import subprocess
import csvutil



def main():

    global data_path
    global data_root
    global data_rel

    global wrk_dir
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    file_list = glob.glob("%s/generated/html/*.html" % (wrk_dir))
    
    
    log_file = "logs/make-htmldb.log"
    msg = "make-htmldb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")

    for in_file in file_list:
        page_id = in_file.split("/")[-1].split(".")[0]
        cn = open(in_file, "r").read();
        doc = {"pageid":page_id, "cn":cn}
        out_file = wrk_dir + "/jsondb/htmldb/%s.json" % (page_id)
        with open(out_file,"w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        msg = "created jsondb/htmldb/%s.json" % (page_id)
        csvutil.write_log_msg(log_file, msg, "a")

    cmd = "chmod -R 775 " + wrk_dir + "/jsondb/htmldb"
    x, y = subprocess.getstatusoutput(cmd)

if __name__ == '__main__':
    main()


