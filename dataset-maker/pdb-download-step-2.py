import os,sys
import string
import glob
import json
import subprocess
import time


###############################
def main():

    count = 0
    out_dir = "downloads/pdb/current/"
    line_list = open("logs/pdb-download-list.txt", "r").read().split("\n")
    for line in line_list[:-1]:
        species, pdb_id = line.strip().split(",")
        out_file = "%s/%s.pdb" % (out_dir, pdb_id)
        cmd = "wget https://files.rcsb.org/download/%s.pdb" % (pdb_id)
        cmd += " -O %s" % (out_file)
        x = subprocess.getoutput(cmd)
        file_size = os.path.getsize(out_file)
        if file_size == 0:
            cmd = "mv %s downloads/pdb/notfound/" % (out_file)
            x = subprocess.getoutput(cmd)
        count += 1
        if count%20 == 0:
            time.sleep(5)
    

         

if __name__ == '__main__':
	main()

