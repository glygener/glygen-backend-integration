import os,sys
import string
import subprocess
from optparse import OptionParser
import glob
import json
import subprocess
import csvutil




def main():

    global jsondb_dir
    global wrk_dir 


    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker"
    jsondb_dir = wrk_dir + "/jsondb/"

    row_count = 0
    batch_size = 1000000
    prev_grp = 1
    out_file = "jsondb/urldb/grp.1.csv"
    FW = open(out_file, "w")
    seen = {}
    for record_type in ["protein", "glycan", "biomarker"]:
        file_list = glob.glob(jsondb_dir + "%sdb/*.json" % (record_type))
        for in_file in file_list:
            cmd = "grep http %s | grep url | sort -u" % (in_file)
            line_list = subprocess.getoutput(cmd).split("\n")
            for line in line_list:
                for url in line.strip().split(" "):
                    idx = url.find("//")
                    if idx != -1:
                        url = url.replace(",", "").replace("\"", "")
                        if url not in seen:
                            domain = url[idx:].split("/")[1]
                            grp = int(row_count/batch_size) + 1
                            if grp != prev_grp:
                                FW.close()
                                out_file = "jsondb/urldb/grp.%s.csv" % (grp)
                                FW = open(out_file, "w")
                            FW.write("%s\n" % (url))
                            seen[url] = True
                            row_count += 1
                            prev_grp = grp
    FW.close()
 






if __name__ == '__main__':
	main()



