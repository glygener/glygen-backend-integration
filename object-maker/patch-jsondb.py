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


    if True:
        file_list = glob.glob("jsondb/glycandb/*.json")
        for in_file in file_list:
            record_id = in_file.split("/")[-1].replace(".json","")
            doc = json.load(open(in_file))
            flag = False
            for obj in doc["publication"]:
                tmp_list = []
                for o in obj["evidence"]:
                    if o["database"] == "" and o["url"] == "" and o["id"] != "":
                        flag = True
                    else:
                        tmp_list.append(o)
                obj["evidence"] = tmp_list
            if flag:
                out_file = "jsondb/patchdb/%s.json" % (record_id)
                print (out_file)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(doc, indent=4)))



    if False:
        file_list = glob.glob("jsondb/proteindb/*.json")
        for in_file in file_list:
            canon = in_file.split("/")[-1].replace(".json","")
            doc = json.load(open(in_file))
            flag = False
            for obj in doc["crossref"]:
                if obj["database"] == "UniCarbKB" and obj["categories"] == []:
                    obj["categories"] = ["Glycoproteomics"]
                    flag = True
            if flag:
                out_file = "jsondb/patchdb/%s.json" % (canon)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(doc, indent=4)))
          
    
      

        
    return



if __name__ == '__main__':
    main()


