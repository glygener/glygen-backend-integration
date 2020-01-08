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


            field_dict = {}
            in_file = "generated/misc/field_names.csv"
            libgly.load_sheet_as_dict(field_dict, in_file, ",", "field_name_current")

            seen = {}
            for doc in dbh["c_bco"].find({}):
                if "bco_id" not in doc:
                    continue
                bco_id = doc["bco_id"]
                if "io_domain" not in doc:
                    continue
                if doc["io_domain"]["output_subdomain"] == []:
                    continue
               
                status_list = []
                if "dataset_categories" in doc["extension_domain"]:
                    for o in doc["extension_domain"]["dataset_categories"]:
                        if o["category_name"] == "status":
                            status_list.append(o["category_value"].lower())

                if "retired" in status_list:
                    continue


                file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"]
                file_ext = file_name.split(".")[-1]
                in_file = data_dir + "/" + file_name
                stat_file_name = ".".join(file_name.split(".")[0:-1]) + ".stat.csv"
                stat_file = data_dir + "/" + stat_file_name
                if os.path.isfile(in_file) == False or os.path.isfile(stat_file) == False:
                    print in_file
                    print stat_file

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

