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
        
        coll = options.coll

        try:
            client = MongoClient('mongodb://localhost:27017',
                username=db_obj["mongodbuser"],
                password=db_obj["mongodbpassword"],
                authSource=db_obj["mongodbname"],
                authMechanism='SCRAM-SHA-1',
                serverSelectionTimeoutMS=10000
            )
            client.server_info()
            mongo_dbh = client[db_obj["mongodbname"]]
            mongo_cl = mongo_dbh[db_obj["collections"][coll]]
            result = mongo_cl.delete_many({})
            out_json = {"taskstatus":1}
        except pymongo.errors.ServerSelectionTimeoutError as err:
            out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
        except pymongo.errors.OperationFailure as err:
            out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}


        print json.dumps(out_json, indent=4)




if __name__ == '__main__':
	main()

