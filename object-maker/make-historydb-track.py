import os,sys
import string
import subprocess
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
from Bio import SeqIO
import gzip

import csvutil




__version__="1.0"
__status__ = "Dev"







def sort_release_list(tmp_list, reversed_flag):

    factor_list = [10000, 1000, 1]
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


###############################
def main():


    current_rel = "x.x.x"
    wrk_dir = "/home/rykahsay/glygen-backend-integration/object-maker"
    data_release_dir = "/data/shared/glygen/releases/data/"

    url = "https://api.glygen.org//misc/verlist/"
    cmd = "curl -s -k %s" % (url)
    x,res = subprocess.getstatusoutput(cmd)
    tmp_list = json.loads(res)
    release_list = sort_release_list(tmp_list, False)
    release_list.append(current_rel)

    his_dict = {}
    for rel in release_list:
        historydb_dir = data_release_dir + "/v-%s/jsondb/historydb/" % (rel)
        if rel == current_rel:
            historydb_dir = wrk_dir + "/jsondb/historydb/"
        file_list = glob.glob(historydb_dir+"/GLY_*.pairs.json")
        for in_file in file_list:
            bco_id = in_file.split("/")[-1].replace(".pairs.json", "")
            doc = json.loads(open(in_file, "r").read())
            r_list = list(doc["history"].keys())
            r_list = sort_release_list(r_list, False)
            row_count_one , row_count_two = 0, 0
            if len(r_list) == 1:
                row_count_two = doc["history"][r_list[0]]["row_count"]
            else:
                row_count_one = doc["history"][r_list[0]]["row_count"]
                row_count_two = doc["history"][r_list[1]]["row_count"]
            for r in r_list:
                tmp_obj = {"row_count_last": row_count_one, "row_count_change":row_count_two - row_count_one}
                for k in doc["history"][r]:
                    if k in ["fields_added", "fields_removed", "ids_added", "ids_removed"]:
                        tmp_obj[k] = len(doc["history"][r][k])
                    else:
                        tmp_obj[k] = doc["history"][r][k]
                if r == rel.replace(".", "_"):
                    if bco_id not in his_dict:
                        his_dict[bco_id] = {}
                    his_dict[bco_id][r] = tmp_obj
    for bco_id in his_dict:
        tmp_obj = {}
        r_list = sort_release_list(list(his_dict[bco_id].keys()), False)
        for r in r_list:
            tmp_obj[r] = his_dict[bco_id][r]
            rel = r.replace("_",".")
            if rel == current_rel:
                out_file = wrk_dir + "/jsondb/historydb/%s.track.json" % (bco_id)
                doc = {"bcoid":bco_id, "doctype":"track", "history":tmp_obj}
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(doc, indent=4)))
            



if __name__ == '__main__':
	main()

