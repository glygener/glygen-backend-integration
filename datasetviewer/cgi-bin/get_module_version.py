#!/usr/bin/python
import os,sys
import string
import commands
from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"




###############################
def main():

    config_json = json.loads(open("conf/config.json", "r").read())
    custom_config_json = json.loads(open("conf/config.custom.json", "r").read())
    db_obj = custom_config_json[config_json["server"]]["dbinfo"]
    path_obj = custom_config_json[config_json["server"]]["pathinfo"]
    root_obj = custom_config_json[config_json["server"]]["rootinfo"]

    client = MongoClient('mongodb://localhost:27017')
    db = client[db_obj["dbname"]]

    out_json = {"version":""}
    for doc in db["c_version"].find({"component":"data"}):
        doc.pop("_id")
        if "version" in doc:
            doc["version"] = doc["version"].replace("v-", "")
            out_json["version"] = "version-%s" % (doc["version"])
        break
    print json.dumps(out_json, indent=4)


if __name__ == '__main__':
	main()

