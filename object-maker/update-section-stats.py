#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import libgly
import csvutil
import section_stats




def main():


    with open("./logfile", "w") as FL:
        FL.write("")

    file_list = glob.glob("jsondb/proteindb/*.json")
    idx = 0
    n = 0
    for in_file in file_list:
        idx += 1
        canon = in_file.split("/")[-1].replace(".json", "")
        doc = json.loads(open(in_file, "r").read())
        flag = False
        if len(doc["glycosylation"]) > 0:
            flag = True 
        if flag:
            out_file = "temp/proteindb/%s.json" % (canon)
            with open("./logfile", "a") as FL:
                FL.write("%s,%s,%s\n" % (idx, len(file_list), out_file))
            section_stats.get_protein_sec_stats(doc, "protein")
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
    
    return



if __name__ == '__main__':
    main()


