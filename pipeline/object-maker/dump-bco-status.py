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
        reviewed_dir = "/data/projects/glygen/generated/datasets/reviewed/"

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
                cond_list = []
                cond_list.append("bco_id" not in doc)
                cond_list.append("extension_domain" not in doc)
                cond_list.append("io_domain" not in doc)
                if True in cond_list:
                    continue
                cond_list.append("dataset_categories" not in doc["extension_domain"])
                cond_list.append("output_subdomain" not in doc["io_domain"])
                if True in cond_list:
                    continue
                cond_list.append(doc["io_domain"]["output_subdomain"] == [])
                if True in cond_list:
                    continue

                bco_id = doc["bco_id"].split("/")[-1]
                file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"].strip()
                for o in doc["extension_domain"]["dataset_categories"]:
                    name = o["category_name"]
                    value = o["category_value"]
                    if name.lower() == "status":
                        file_path = reviewed_dir + file_name
                        flag = os.path.isfile(file_path)
                        print bco_id,value.lower(), file_name, flag


        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

