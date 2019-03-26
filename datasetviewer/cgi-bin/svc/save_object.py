import os,sys
import string
import commands
from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import pymongo
from pymongo import MongoClient
import traceback


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
       
        out_json = {}
	try:
		client = MongoClient('mongodb://localhost:27017')
		db = client[db_obj["dbname"]]
        	coll = "c_metadata"
        	obj = json.loads(options.injson)
        	obj_id = int(obj["objid"])
        	obj.pop("objid")
		if obj == {}:
			out_json = {"status":-1}
			out_json["errormsg"] = "Cannot save empty object"
			print json.dumps(out_json, indent=4)
			sys.exit()
		elif obj_id == 0:
            		obj["_id"]  = get_next_sequence_value(db, "metadataid")
            		result = db[coll].insert_one(obj)
            		obj_id = obj["_id"]
        	else:
            		result = db[coll].replace_one({"_id":obj_id}, obj, upsert=False)

		query_obj = {"_id":obj_id}
            	out_json = db[coll].find_one(query_obj)
		out_json["status"] = 1
                out_json["errormsg"] = ""
	except Exception, e:
		out_json = {"status":-1}
                out_json["errormsg"] =  traceback.format_exc().split("\n")[-2].split(":")[0]
		#print traceback.format_exc()

        print json.dumps(out_json, indent=4)


if __name__ == '__main__':
	main()

