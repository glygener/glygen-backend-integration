import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


def reset_sequence_value(db, sequence_name):
    seq_doc = db.c_counters.find_and_modify(
        query={'_id': sequence_name},
        update={'$set': {'sequence_value': 0}},
        upsert=True,
        new=True
    )
    return


def get_next_sequence_value(db, sequence_name):

    seq_doc = db.c_counters.find_and_modify(
        query={'_id': sequence_name}, 
        update={'$inc': {'sequence_value': 1}}, 
        upsert=True, 
        new=True
    )
    return int(seq_doc["sequence_value"])




###############################
def main():


        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-i","--infile",action="store",dest="infile",help="In file")
        parser.add_option("-r","--recordid",action="store",dest="recordid",help="Record ID")



        (options,args) = parser.parse_args()
        for key in ([options.infile, options.recordid]):
            if not (key):
                parser.print_help()
                sys.exit(0)
                                                                                                 

        in_file = options.infile
        record_id = options.recordid

        obj = json.loads(open(in_file, "r").read())
        if record_id in obj:
            print json.dumps(obj[record_id], indent=4)
        else:
            print record_id, "not found in json db"



if __name__ == '__main__':
	main()

