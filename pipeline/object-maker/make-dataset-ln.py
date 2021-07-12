import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
from Bio import SeqIO
import csv



__version__="1.0"
__status__ = "Dev"



###############################
def main():

        config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())
        db_obj = config_obj["dev"]["dbinfo"]
        data_dir = "reviewed/"
        coll = "c_bco"
        #data_dir = "/data/shared/glygen/releases/data/v-1.5.36/reviewed/"
        #coll = "c_bco_v-1.5.36"

        os.chdir(data_dir)


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

            
            seen = {}
            for doc in dbh[coll].find({}):
                if "bco_id" not in doc:
                    continue
                bco_id = doc["bco_id"]
                bco_num = bco_id.split("/")[-1]
                if "io_domain" not in doc:
                    continue
                if doc["io_domain"]["output_subdomain"] == []:
                    continue
               
                retired_flag = "active"
                if "dataset_categories" in doc["extension_domain"]:
                    for o in doc["extension_domain"]["dataset_categories"]:
                        if o == None:
                            continue
                        cat_name = o["category_name"].replace(" ", "_")
                        cat_value = o["category_value"].strip()
                        if cat_name in ["status", "dataset_status"] and cat_value == "retired":
                            retired_flag = "retired"
                file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"].strip()
                if os.path.isfile(file_name):
                    file_ext = file_name.split(".")[-1]
                    cmd = "ln -s %s %s.%s" % (file_name, bco_num, file_ext)
                    x = commands.getoutput(cmd)

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

