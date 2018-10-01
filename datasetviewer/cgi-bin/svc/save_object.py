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

def get_next_sequence_value(db, sequence_name):
        
    seq_doc = db.c_counters.find_and_modify(
        query={'_id': sequence_name}, 
        update={'$inc': {'sequence_value': 1}}, 
        upsert=True, 
        new=True
    )
    return int(seq_doc["sequence_value"])




###############################
def main():


        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-j","--injson",action="store",dest="injson",help="Object ID")
        
        (options,args) = parser.parse_args()
        for key in ([options.injson]):
            if not (key):
                parser.print_help()
                sys.exit(0)



        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
	client = MongoClient('mongodb://localhost:27017')
	db = client[db_obj["dbname"]]
        coll = "c_metadata"
        
        obj = json.loads(options.injson)
        obj_id = int(obj["objid"])

        if obj_id == 0:
            obj.pop("objid")
            obj["_id"]  = get_next_sequence_value(db, "metadataid")
            result = db[coll].insert_one(obj)
            obj_id = obj["_id"]
        else:
            result = db[coll].update_one({"_id":obj_id}, {'$set': obj}, upsert=False)
   
        query_obj = {"_id":obj_id}
        #doc = db[coll].find_one(query_obj)
        #doc.pop("_id")
        #doc.pop("object_id")

        doc = {"objid":obj_id}
        print json.dumps(doc, indent=4)


if __name__ == '__main__':
	main()

