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

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
	client = MongoClient('mongodb://localhost:27017')
	db = client[db_obj["dbname"]]
        coll = "c_metadata"
        query_obj = {}

        out_json = {"datasets":[]}
        for doc in db[coll].find(query_obj).sort("_id", 1):
            doc = json.loads(json_util.dumps(doc))
            out_json["datasets"].append(doc)
        
        print json.dumps(out_json, indent=4)
                        


if __name__ == '__main__':
	main()

