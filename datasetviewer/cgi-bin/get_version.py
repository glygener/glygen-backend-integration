#!/usr/bin/python
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
import cgi


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
    
    print "Content-Type: application/json"
    print   

    out_json = {}
    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        version_one, version_two = config_json["moduleversion"],config_json["datarelease"]
        versions = "Viewer-%s | Data-%s" % (version_one, version_two)
        out_json = {"moduleversion":versions, "taskstatus":1};
    except Exception, e:
        out_json = {"taskstatus":0, "errormsg":"query failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

