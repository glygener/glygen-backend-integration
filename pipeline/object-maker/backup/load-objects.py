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


def reset_sequence_value(mongo_cl_counters, sequence_name):
    seq_doc = mongo_cl_counters.find_and_modify(
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
        
        (options,args) = parser.parse_args()
        for key in ([options.coll]):
            if not (key):
                parser.print_help()
                sys.exit(0)
                                                                                                 
        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
        coll = options.coll
        
        try:
            client = MongoClient('mongodb://localhost:27017',
                username=db_obj["mongodbuser"],
                password=db_obj["mongodbpassword"],
                authSource=db_obj["mongodbname"],
                authMechanism='SCRAM-SHA-1',
                serverSelectionTimeoutMS=10000
            )
            client.server_info()
            mongo_dbh = client[db_obj["mongodbname"]]
            mongo_cl_counters = mongo_dbh[db_obj["collections"]["counters"]]
        
            mongo_cl = mongo_dbh[db_obj["collections"][coll]]

            nrecords = 0
            if coll == "glycan":
                primary_id = "glycanid"
                json_db_file = path_obj["jsondbpath"] + "glycandb.json"
                objs = json.loads(open(json_db_file, "r").read()) 
                reset_sequence_value(mongo_cl_counters, primary_id)
                for key in objs:
                    objs[key]["_id"]  = get_next_sequence_value(mongo_cl_counters, primary_id)
                    result = mongo_cl.insert_one(objs[key])
                    nrecords += 1
            if coll == "protein":
                primary_id = "proteinid"
                json_db_file = path_obj["jsondbpath"] + "proteindb.json" 
                objs = json.loads(open(json_db_file, "r").read())
                reset_sequence_value(mongo_cl_counters, primary_id) 
                for key in objs:
                    objs[key]["_id"]  = get_next_sequence_value(mongo_cl_counters, primary_id)
                    result = mongo_cl.insert_one(objs[key])
                    nrecords += 1
            elif coll == "searchinit":
                json_db_file = path_obj["jsondbpath"] + "searchinitdb.json"
                objs = json.loads(open(json_db_file, "r").read())
                result = mongo_cl.insert_one(objs)
                nrecords += 1
            out_json = {"taskstatus":1, "objects_loaded":nrecords}
        except pymongo.errors.ServerSelectionTimeoutError as err:
            out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
        except pymongo.errors.OperationFailure as err:
            out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}


        print json.dumps(out_json, indent=4)









if __name__ == '__main__':
	main()

