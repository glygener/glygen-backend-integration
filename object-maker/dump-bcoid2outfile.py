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



def load_bcos_from_editor(bco_dict, filename2bcoid, bcoid2filename): 
    headers = {
        "Authorization":"Token 5467202e22e657dda9c19edb9a25936125163373",
        "Content-Type":"application/json"
    }
    req_obj = {
        "POST_api_objects_search":[
            {"type":"bco_id","search":"GLY_"}
        ]    
    }
    
    cmd =  'curl --request POST https://biocomputeobject.org/api/objects/search/'
    cmd += ' -H ' + '"Content-type:application/json"'
    cmd += ' -H ' + '"Authorization:Token 5467202e22e657dda9c19edb9a25936125163373"'
    cmd += ' -d ' + '\'{"POST_api_objects_search":[{"type":"bco_id","search":"GLY_"}]}\''
    cmd += ' -o logs/bco.txt'
    x = subprocess.getoutput(cmd)
    
    x = subprocess.getoutput("cat logs/bco.txt")
    bco_obj_list = json.loads(x)
    for bco_doc in bco_obj_list:
        bco_id = bco_doc["object_id"].split("/")[-2]

        bco_ver = bco_doc["object_id"].split("/")[-1]
        if bco_ver == "DRAFT" and "contents" in bco_doc:
            status_list = []            
            if "extension_domain" in bco_doc["contents"]:
                for obj in bco_doc["contents"]["extension_domain"]:
                    if "dataset_extension" in obj:
                        o = obj["dataset_extension"]
                        if "dataset_categories" in o:
                            for oo in o["dataset_categories"]:
                                cat_name, cat_value = oo["category_name"].lower(), oo["category_value"].lower()
                                {'category_value': 'Reviewed', 'category_name': 'status'}
                                if cat_name == "status":
                                    status_list.append(cat_value)
            if "retired" in status_list:
                continue
            bco_dict[bco_id] = bco_doc["contents"]
            if "io_domain" in bco_doc["contents"]:
                for obj in bco_doc["contents"]["io_domain"]["output_subdomain"]:
                    if "uri" in obj:
                        if "filename" in obj["uri"]:
                            file_name = obj["uri"]["filename"]
                            if file_name.strip() != "" and file_name.find(".stat.csv") == -1:
                                filename2bcoid[file_name] = bco_id
                                if bco_id not in bcoid2filename:
                                    bcoid2filename[bco_id] = {}
                                bcoid2filename[bco_id][file_name] = True
        
    return 







def main():


    bco_dict, fname2bcoid, bcoid2fname = {}, {}, {}
    load_bcos_from_editor(bco_dict, fname2bcoid, bcoid2fname)

    for bco_id in bcoid2fname:
        for file_name in bcoid2fname[bco_id]:
            print (bco_id, file_name)





if __name__ == '__main__':
    main()


