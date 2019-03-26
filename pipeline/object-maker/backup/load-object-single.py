import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


def reset_sequence_value(db, sequence_name):
    seq_doc = db.c_counters.find_and_modify(
        query={'_id': sequence_name},
        update={'$set': {'sequence_value': 0}},
        upsert=True,
        new=True
    )
    return


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
        parser.add_option("-c","--coll",action="store",dest="coll",help="Collection name")
        parser.add_option("-i","--infile",action="store",dest="infile",help="In file")



        (options,args) = parser.parse_args()
        for key in ([options.coll, options.infile]):
            if not (key):
                parser.print_help()
                sys.exit(0)
                                                                                                 

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
	client = MongoClient('mongodb://localhost:27017')
	db = client[db_obj["dbname"]]
        coll = options.coll
        in_file = options.infile
       
        obj = json.loads(open(in_file, "r").read())
       
        primary_dict = {
            "c_glycan":"glycanid",
            "c_protein":"proteinid",
            "c_searchinit":"searchinitid",
            "c_schema":"schemaid"
        }

        if coll in primary_dict:
            obj["_id"]  = get_next_sequence_value(db, primary_dict[coll])
        
        #print json.dumps(obj, indent=4)
        result = db[coll].insert_one(obj)	


if __name__ == '__main__':
	main()

