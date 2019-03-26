import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient

sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"



###############################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        data_dir = "reviewed/"

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


            field_dict = {}
            in_file = data_dir + "/field_names.csv"
            libgly.load_sheet_as_dict(field_dict, in_file, ",", "field_name_current")

            seen = {}
            for doc in dbh["c_bco"].find({}):
                if "bco_id" not in doc:
                    continue
                query_obj = {"bco_id":doc["bco_id"]}
                update_obj = {
                    "error_domain":{
                        "emperical_error": {},
                        "algorithmic_error":{}
                    }
                }
                result = dbh["c_bco"].update_one(query_obj, {'$set': update_obj}, upsert=True)
                print "updated %s" % (doc["bco_id"])

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

