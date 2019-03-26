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

import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"




###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-o","--objid",action="store",dest="objid",help="Object ID")
    parser.add_option("-v","--objver",action="store",dest="objver",help="Object Version")


    (options,args) = parser.parse_args()
    for key in ([options.objid]):
        if not (key):
            parser.print_help()
            sys.exit(0)


    global config_json
    global db_obj
    global client
    global path_obj


    out_json = {}

    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        db_obj = config_json[config_json["server"]]["dbinfo"]
        path_obj = config_json[config_json["server"]]["pathinfo"]

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
        obj_id = "http://data.glygen.org/%s" % (options.objid)
        
        current_ver = open(path_obj["htmlpath"] + "/release-notes.txt", "r").read().strip().split(" ")[-1]
        obj_ver = options.objver if options.objver != None else current_ver

        query_obj = {"bco_id":obj_id}
        bco_collection = "c_bco_" + obj_ver
        doc = dbh[bco_collection].find_one(query_obj)
        ordr_dict = json.loads(open("conf/field_order.json").read())
        doc = util.order_json_obj(doc, ordr_dict)

        print "<pre>"
        
        if doc == None:
            print "\n\n\tReadme file does not exist for object %s!" % (obj_id)
        else:
            doc.pop("_id")
            print json.dumps(doc, indent=4)
        print "</pre>"
        sys.exit()
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

