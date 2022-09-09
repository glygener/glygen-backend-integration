import os,sys
import string
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"




def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-c","--coll",action="store",dest="coll",help="Collection")
        parser.add_option("-v","--ver",action="store",dest="ver",help="Release version")

        (options,args) = parser.parse_args()
        for key in ([options.coll, options.ver]):
            if not (key):
                parser.print_help()
                sys.exit(0)

        ver = options.ver
        coll = options.coll
        db_obj = {
            "mongodbuser":"***",
            "mongodbpassword":"***",
            "mongodbname":"tmpdb"
        }


        db_path = "/data/shared/glygen/releases/data/v-%s/jsondb/" % (ver)

        file_list = glob.glob(db_path + "%sdb/*.json" % (coll.replace("c_","")))
        
        #print (file_list)
        #exit()


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
            nrecords = 0
            for in_file in file_list:
                doc = json.loads(open(in_file, "r").read())
                if "_id" in doc:
                    doc.pop("_id")
                c = coll
                if coll in ["c_bco", "c_extract", "c_history"]:
                    c = "%s_v-%s" % (coll, ver)
                result = dbh[c].insert_one(doc)     
                nrecords += 1
            print (" ... loaded %s objects to %s" % (nrecords, c))
        except pymongo.errors.ServerSelectionTimeoutError as err:
            print (err)
        except pymongo.errors.OperationFailure as err:
            print (err)



if __name__ == '__main__':
        main()
