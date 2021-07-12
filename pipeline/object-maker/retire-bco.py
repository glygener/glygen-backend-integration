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


        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i","--bcoid",action="store",dest="bcoid",help="")

        (options,args) = parser.parse_args()
        for file in ([options.bcoid]):
            if not (file):
                parser.print_help()
                sys.exit(0)
        

        bco_number = int(options.bcoid)


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
                bco_idx = bco_id.split("_")[-1]
                if bco_idx.find("id") != -1:
                    continue
                bco_idx = int(bco_idx)
                if bco_idx != bco_number:
                    continue
                for o in doc["extension_domain"]["dataset_categories"]:
                    name = o["category_name"]
                    value = o["category_value"]
                    if name == "status":
                        o["category_value"] = "Retired"
                update_obj = {"extension_domain":doc["extension_domain"]}
                query_obj = {"bco_id":bco_id}
                r = dbh["c_bco"].update_one(query_obj, {'$set': update_obj}, upsert=True)
                #print json.dumps(update_obj, indent=4)
                #print json.dumps(query_obj, indent=4)

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

