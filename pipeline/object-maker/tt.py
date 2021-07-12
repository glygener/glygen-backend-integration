import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import datetime
import pymongo
from pymongo import MongoClient
from Bio import SeqIO

sys.path.append('../../glytools/')
import libgly


__version__="1.0"
__status__ = "Dev"



###############################
def main():

    obj_list = {}
    for in_file in glob.glob("jsondb/pathdb/*.json"):
        record_type = in_file.split("/")[-1].split(".")[0]
        doc = json.loads(open(in_file, "r").read())
        path, data_type = doc["path"].split(" ")
        doc["path"] = path
        doc["type"] = data_type

        if record_type not in obj_list:
            obj_list[record_type] = []
        obj_list[record_type].append(doc)
    
    print json.dumps(obj_list, indent=4)


if __name__ == '__main__':
	main()

