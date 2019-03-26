import os,sys
import string
import commands
import csv
import traceback

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import util

import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"




###############################
def main():

    global config_json
    global db_obj
    global client
    global path_obj


    out_json = {}

    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        db_obj = config_json[config_json["server"]]["dbinfo"]
        path_obj = config_json[config_json["server"]]["pathinfo"]

    except Exception, e:
        out_json = {"taskstatus":0, "errormsg":"Loading config failed!"}

    try:
        client = MongoClient('mongodb://localhost:27017',
            username=db_obj["mongodbuser"],
            password=db_obj["mongodbpassword"],
            authSource=db_obj["mongodbname"],
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        dbh = client[db_obj["mongodbname"]]
        current_ver = open(path_obj["htmlpath"] + "/release-notes.txt", "r").read().strip().split(" ")[-1]

        bco_collection = "c_bco"
        for doc in dbh[bco_collection].find({}):
            if "bco_id" in doc:
                query_obj = {"bco_id":doc["bco_id"]}
                bco_id = doc["bco_id"].replace("GLYDS", "DSBCO_")
                update_obj = {"bco_id":bco_id}
                res = dbh[bco_collection].update_one(query_obj, {'$set': update_obj}, upsert=True)

        sys.exit()
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

