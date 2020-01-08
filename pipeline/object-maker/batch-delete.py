import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
import datetime
import util
import hashlib

sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"


def batch_delete(dbh, offset):


    #Delete object with bco_id > offset
    for doc in dbh["c_bco"].find({}):
        if "bco_id" in doc:
            bco_id = doc["bco_id"]
            bco_id_int = int(bco_id.split("_")[-1])
            if bco_id_int  > offset:
                dbh["c_bco"].delete_one({"bco_id":bco_id})
                print "Deleted ", bco_id_int
    
    #reset bco ID counter
    query_obj = {'_id': "bcoid"}
    update_obj = {'$set': {'sequence_value': offset}}
    seq_doc = dbh["c_counters"].find_and_modify(query=query_obj,update=update_obj,upsert=True,new=True)
    return



###############################
def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-o","--offset",action="store",dest="offset",help="Offset")

        (options,args) = parser.parse_args()
        for key in ([options.offset]):
            if not (key):
                parser.print_help()
                sys.exit(0)

        offset = int(options.offset)

        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        data_dir = "unreviewed/"

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

            #Batch delete
            batch_delete(dbh, offset)


        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

