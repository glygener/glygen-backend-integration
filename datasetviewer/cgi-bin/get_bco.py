#!/usr/bin/python
import os,sys
import string
import commands
import csv
import traceback

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import util
import cgi


import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"




###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    msg = "Input JSON text"
    parser.add_option("-j","--injson",action="store",dest="injson",help=msg)


    form_dict = cgi.FieldStorage()
    (options,args) = parser.parse_args()
    local_flag = False
    in_json = {}
    if len(form_dict.keys()) > 0:
        in_json = json.loads(form_dict["injson"].value) if "injson" in form_dict else {}
    else:
        local_flag = True
        for key in ([options.injson]):
            if not (key):
                parser.print_help()
                sys.exit(0)
        in_json = json.loads(options.injson)


    global config_json
    global db_obj
    global client
    global path_obj
    
    print "Content-Type: application/json"
    print   

    #print json.dumps(in_json, indent=4)
    #sys.exit()


    out_json = {}

    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        custom_config_json = json.loads(open("conf/config.custom.json", "r").read())
        db_obj = custom_config_json[config_json["server"]]["dbinfo"]
        path_obj = custom_config_json[config_json["server"]]["pathinfo"]
        root_obj = custom_config_json[config_json["server"]]["rootinfo"]
    except Exception, e:
        out_json = {"taskstatus":0, "errormsg":"Loading config failed!"}

    try:
        client = MongoClient('mongodb://localhost:27017',
            username=db_obj["mongodbuser"],
            password=db_obj["mongodbpassword"],
            authSource=db_obj["mongodbname"],
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        dbh = client[db_obj["mongodbname"]]

        in_json["objid"] = in_json["objid"].replace("/", "")

        obj_id = db_obj["bcourl"] % (in_json["objid"].split("_")[1]) 

        current_ver = open(path_obj["htmlpath"] + "/release-notes.txt", "r").read().strip().split(" ")[-1]
        obj_ver = "v-" + config_json["datarelease"]
        if "objver" in in_json:
            obj_ver = in_json["objver"] if in_json["objver"] != None else obj_ver


        query_obj = {"bco_id":obj_id}
        bco_collection = "c_bco_" + obj_ver
        doc = dbh[bco_collection].find_one(query_obj)
        

        if doc == None:
            print "\n\n\tReadme file does not exist for object %s!" % (obj_id)
        else:
            ordr_dict = json.loads(open("conf/field_order.json").read())
            doc = util.order_json_obj(doc, ordr_dict)
            if "_id" in doc:
                doc.pop("_id")
            print json.dumps(doc, indent=4)
        sys.exit()
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

