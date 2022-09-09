import os,sys
import string
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"



###############################
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

        coll = options.coll
        ver = options.ver

        if coll in ["c_bco", "c_history", "c_extract"]:
            coll = "%s_v-%s" % (coll, ver)
       
        db_obj = {
            "mongodbuser":"***",
            "mongodbpassword":"***",
            "mongodbname":"tmpdb"
        }

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
            result = dbh[coll].delete_many({})
            #print ("Cleaned collection %s" % (coll))
       
        except pymongo.errors.ServerSelectionTimeoutError as err:
            print (err)
        except pymongo.errors.OperationFailure as err:
            print (err)



if __name__ == '__main__':
	main()
