import os,sys
import string
import glob
import json
import subprocess
import time


###############################
def main():


    count = 0
    seen = {}
    out_dir = "downloads/pdb/current/"
    file_list = glob.glob("unreviewed/*_protein_pdb_map.csv")
    for in_file in file_list:
        species = in_file.split("/")[-1].split("_")[0]
        line_list = open(in_file, "r").read().split("\n")
        for line in line_list[1:-1]:
            row = line.split(",")
            method = row[-3].replace("\"","")
            if method == "AlphaFold":
                continue
            pdb_id = row[2].replace("\"","").lower() 
            selected = row[-1].replace("\"","")
            if pdb_id not in seen and selected == "True":
                out_file = "%s/%s.pdb" % (out_dir, pdb_id)
                if os.path.isfile(out_file) == False:
                    continue
                file_size = os.path.getsize(out_file)
                if file_size == 0:
                    print (pdb_id)
                seen[pdb_id] = True
 

if __name__ == '__main__':
	main()

