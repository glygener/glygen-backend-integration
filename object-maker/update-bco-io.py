#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import requests
import subprocess
import pymongo
from optparse import OptionParser
import libgly
from Bio import SeqIO


__version__="1.0"
__status__ = "Dev"



def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-i","--bcoid",action="store",dest="bcoid",help="BCO ID")
    parser.add_option("-o","--olduri",action="store",dest="olduri",help="Existing file URI")
    parser.add_option("-n","--newuri",action="store",dest="newuri",help="New file URI")
    parser.add_option("-f","--filename",action="store",dest="filename",help="New file name")


    (options,args) = parser.parse_args()
    for key in ([options.bcoid]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    bco_id = options.bcoid
    old_uri = options.olduri
    new_uri = options.newuri
    new_filename = options.filename

    db_obj = json.loads(open("../conf/db.json", "r").read())
    try:
        client = pymongo.MongoClient('mongodb://localhost:27017',
            username=db_obj["mongodbuser"],
            password=db_obj["mongodbpassword"],
            authSource=db_obj["mongodbname"],
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        dbh = client[db_obj["mongodbname"]]
        coll = "c_bco"
        qry_obj = {"bco_id":{"$regex":bco_id}}
        doc = dbh[coll].find_one(qry_obj)
        bco_url = doc["bco_id"]

        if "io_domain" in doc:
            if "output_subdomain" in doc["io_domain"]:
                out_file_list = [] 
                for obj in doc["io_domain"]["output_subdomain"]:
                    file_name = obj["uri"]["filename"].strip()
                    ignore = False
                    for k in [".stat.csv", "GLY_0"]:
                        if file_name.find(k) != -1:
                            ignore = True
                    if ignore == True or file_name == "":
                        continue
                    out_file_list.append(file_name)
                update_flag = False
                for obj in doc["io_domain"]["input_subdomain"]:
                    file_name = obj["uri"]["filename"].strip()
                    if old_uri == None:
                        print (file_name, obj["uri"]["uri"])
                    elif obj["uri"]["uri"] == old_uri:
                        if new_uri != None:
                            obj["uri"]["uri"] = new_uri
                            update_flag = True
                        if new_filename != None:
                            obj["uri"]["filename"] = new_filename
                            update_flag = True
                
                if update_flag == True:
                    update_obj = {"io_domain.input_subdomain":doc["io_domain"]["input_subdomain"]}
                    query_obj = {"bco_id":bco_url}
                    r = dbh[coll].update_one(query_obj, {'$set': update_obj}, upsert=True)


    except pymongo.errors.ServerSelectionTimeoutError as err:
        print ({"error_list":[{"error_code": "open-connection-failed"}]})
    except pymongo.errors.OperationFailure as err:
        print ({"error_list":[{"error_code": "mongodb-auth-failed"}]})
    


if __name__ == '__main__':
    main()



