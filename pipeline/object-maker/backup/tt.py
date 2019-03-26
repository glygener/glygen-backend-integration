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


##################################
def parse_prop_tree(prop_json, obj, parent):
    if parent not in prop_json:
        prop_json[parent] = {"children":[]}

    if "properties" in obj:
        for p in obj["properties"]:
            o = obj["properties"][p]
            tmp_obj = {"name":p, "parent":parent}
            if "type" in o:
                tmp_obj["type"] = o["type"]
            prop_json[parent]["children"].append(p)
            parse_prop_tree(prop_json, obj["properties"][p], p)
        




#######################
def parse_obj(obj, parent, child_list):

    if parent not in child_list:
        child_list[parent] = []

    if type(obj) is dict:
        key_list = obj.keys()
        for k1 in key_list:
            field = parent + "." + k1
            val_type = ""
            val_type = "dict" if type(obj[k1]) is dict else val_type
            val_type = "list" if type(obj[k1]) is list else val_type
            print field, val_type
            child_list[parent].append(k1)
            if type(obj[k1]) in [dict, list]:
                parse_obj(obj[k1], parent + ":" + k1, child_list)
    elif type(obj) is list:
        for k1 in xrange(0, len(obj)):
            field = parent  
            val_type = ""
            val_type = "dict" if type(obj[k1]) is dict else val_type
            val_type = "list" if type(obj[k1]) is list else val_type
            print field, val_type
            child_list[parent].append(k1)
            if type(obj[k1]) in [dict, list]:
                parse_obj(obj[k1], parent + ":" + k1, child_list)

    return





###############################
def main():


        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-c","--coll",action="store",dest="coll",help="Collection name")



        (options,args) = parser.parse_args()
        for key in ([options.coll]):
            if not (key):
                parser.print_help()
                sys.exit(0)
                                                                                                 

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
	client = MongoClient('mongodb://localhost:27017')
	db = client[db_obj["dbname"]]
        coll = options.coll
       
        doc = db[coll].find_one({"_id":4})
        prop_json = {}
        #lineage = parse_prop_tree(prop_json, doc, "root")
        child_list = {}
        parse_obj(doc, "root", child_list)
        #for k in child_list:
        #    print k

        print json.dumps(prop_json, indent=4)


if __name__ == '__main__':
	main()

