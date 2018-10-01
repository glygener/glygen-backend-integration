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


        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-o","--objid",action="store",dest="objid",help="Object ID")
        
        (options,args) = parser.parse_args()
        for key in ([options.objid]):
            if not (key):
                parser.print_help()
                sys.exit(0)



        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
	client = MongoClient('mongodb://localhost:27017')
	db = client[db_obj["dbname"]]
        coll = "c_metadata"

        obj_id = int(options.objid.replace("GLYDS", ""))
       
        query_obj = {"_id":obj_id}
        doc = db[coll].find_one(query_obj)

        print "<pre>"
        if doc == None:
            print "\n\n\tReadme file does not exist for object %s!" % (obj_id)
        else:
            readme_file = path_obj["htmlpath"] + "/datasets/reviewed/" + doc["readmefilename"]
            if os.path.exists(readme_file):
                with open(readme_file, "r") as FR:
                    for line in FR:
                        print line[:-1]
            else:
                print "\n\n\tReadme file does not exist for object %s!" % (obj_id)


        print "</pre>"

if __name__ == '__main__':
	main()

