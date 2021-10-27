#!/usr/bin/python
import os,sys
import string
import commands
from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"




###############################
def main():


        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-o","--objid",action="store",dest="objid",help="Object ID")
        
        (options,args) = parser.parse_args()
        for key in ([options.objid]):
            if not (key):
                parser.print_help()
                sys.exit(0)


        config_json = json.loads(open("conf/config.json", "r").read())
        custom_config_json = json.loads(open("conf/config.custom.json", "r").read())
        db_obj = custom_config_json[config_json["server"]]["dbinfo"]
        path_obj = custom_config_json[config_json["server"]]["pathinfo"]
        root_obj = custom_config_json[config_json["server"]]["rootinfo"]

	client = MongoClient('mongodb://localhost:27017')
        db = client[db_obj["dbname"]]
        coll = "c_metadata"
        obj_id = int(options.objid)
       
        if obj_id == 0:
            doc = json.loads(open("conf/template.json", "r").read())
        else:
            query_obj = {"_id":obj_id}
            doc = db[coll].find_one(query_obj)
            doc.pop("_id")
            if "object_id" in doc:
                doc.pop("object_id")
        print json.dumps(doc, indent=4)
                        


if __name__ == '__main__':
	main()

