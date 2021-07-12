import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
from Bio import SeqIO

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

            cat_dict = {}
            for doc in dbh["c_bco"].find({}):
                if "bco_id" not in doc:
                    continue
                if "extension_domain" not in doc:
                    continue
                if "dataset_categories" not in doc["extension_domain"]:
                    continue
                bco_id = doc["bco_id"]
                for o in doc["extension_domain"]["dataset_categories"]:
                    name = o["category_name"]
                    value = o["category_value"]
                    if name == "tags":
                        continue
                    if name not in cat_dict:
                        cat_dict[name] = []
                    if value not in cat_dict[name]:
                        cat_dict[name].append(value)
            print json.dumps(cat_dict, indent=4)

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

