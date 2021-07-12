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



###############################
def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-d","--dbname",action="store",dest="dbname",help="Database name")

        (options,args) = parser.parse_args()
        for key in ([options.dbname]):
            if not (key):
                parser.print_help()
                sys.exit(0)

        db_name = options.dbname
        server = "dev"
        
        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[server]["dbinfo"]


        try:
            client = pymongo.MongoClient('mongodb://localhost:27017',
                username=db_obj["mongodbuser"],
                password=db_obj["mongodbpassword"],
                authSource=db_name,
                authMechanism='SCRAM-SHA-1',
                serverSelectionTimeoutMS=10000
            )
            client.server_info()
            dbh = client[db_name]
            coll = "c_bco_v-1.5.36" 
            #coll = "c_bco"
            for doc in dbh[coll].find({}):
                if "bco_id" not in doc:
                    continue
                if doc["bco_id"].find("DSBCO_") == -1:
                    continue
                q_obj = {"bco_id":doc["bco_id"]}
                newbco_id = doc["bco_id"].replace("DSBCO_", "GLY_")
                update_obj = {"bco_id":newbco_id}
                result = dbh[coll].update_one(q_obj, {'$set': update_obj}, upsert=True)
                print doc["bco_id"], update_obj
                

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

