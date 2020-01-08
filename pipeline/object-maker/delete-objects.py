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



###############################
def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-d","--dbname",action="store",dest="dbname",help="Database name")
        parser.add_option("-c","--coll",action="store",dest="coll",help="Collection name")

        (options,args) = parser.parse_args()
        for key in ([options.dbname, options.coll]):
            if not (key):
                parser.print_help()
                sys.exit(0)

        db_name = options.dbname
        server = "beta" if db_name == "glydb_beta" else "dev"
        
        coll = options.coll
        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[server]["dbinfo"]


        try:
            client = pymongo.MongoClient('mongodb://localhost:27017',
                username=db_obj["mongodbuser"],
                password=db_obj["mongodbpassword"],
                authSource=db_name,
                authMechanism='SCRAM-SHA-1',
                serverSelectionTimeoutMS=10000
            )
            client.server_info()
            dbh = client[db_name]
            result = dbh[coll].delete_many({})
            print "Done!"
        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

