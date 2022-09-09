import os,sys
import datetime,time
import pytz
import json

from optparse import OptionParser
import pymongo

__version__="1.0"
__status__ = "Dev"


###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-d","--dbname",action="store",dest="dbname",help="Database")
    parser.add_option("-v","--compversion",action="store",dest="compversion",help="version")
    parser.add_option("-c","--component",action="store",dest="component",help="component")

    (options,args) = parser.parse_args()
    for key in ([options.dbname, options.compversion, options.component]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    db_name = options.dbname   
    version = options.compversion
    component = options.component
    ts = datetime.datetime.now(pytz.timezone('US/Eastern')).strftime('%Y-%m-%d %H:%M:%S %Z%z')
   

    user_info = {
        "glydb":{"mongodbuser":"***" ,"mongodbpassword":"***"},
        "glydb_beta":{"mongodbuser":"***","mongodbpassword":"***"}
    } 


    out_json = {} 
    try:
        client = pymongo.MongoClient('mongodb://localhost:27017',
            username=user_info[db_name]["mongodbuser"],
            password=user_info[db_name]["mongodbpassword"],
            authSource=db_name,
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        mongo_dbh = client[options.dbname]
        mongo_cl = mongo_dbh["c_version"]
        ver_obj = { "component":component, "version":version, "release_date":ts}
        result = mongo_cl.update_one({"component":component}, {'$set': ver_obj}, upsert=True)
        doc = mongo_cl.find_one({"component":component})
        c, v, d = doc["component"], doc["version"], doc["release_date"]
        print (" ... %s version updated to %s on %s"  % (c, v, d))
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
        print (err)
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}
        print (err)




            

if __name__ == '__main__':
	main()

