import os,sys
import string
import glob
import json
import subprocess
import requests
import time

###############################
def main():

    out_dir = "downloads/pdb/current/"


    count = 0
    line_list = open("logs/alphafold-download-list.txt", "r").read().split("\n")
    for line in line_list[:-1]:
        species, pdb_id = line.strip().split(",")
        ac = pdb_id.split("-")[1]
        url = "https://alphafold.ebi.ac.uk/api/prediction/%s" % (ac)
        res = requests.get(url)
        if res.content.strip() != "":
            res_json = json.loads(res.content)
            for obj in res_json:
                pdb_url = obj["pdbUrl"] if "pdbUrl" in obj else ""
                out_file = "%s/%s.pdb" % (out_dir, pdb_id)
                cmd = "wget %s -O %s" % (pdb_url, out_file)
                x = subprocess.getoutput(cmd)
                file_size = os.path.getsize(out_file)
                if file_size == 0:
                    cmd = "mv %s downloads/pdb/notfound/" % (out_file)
                    x = subprocess.getoutput(cmd)
                count += 1
                if count%10 == 0:
                    time.sleep(5)


if __name__ == '__main__':
	main()

