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
    
    config_obj = json.loads(open("conf/config.json", "r").read())
    db_obj = config_obj["dev"]["dbinfo"]
    db_path = "jsondb/"
    db_name = "glydb"


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
        query_dict = {
            "c_glycan":{"glytoucan_ac":"G17689DH"},
            "c_protein":{"uniprot_canonical_ac":"P14210-1"}
        }
        for coll in query_dict:
            query_obj = query_dict[coll]
            doc = dbh[coll].find_one(query_obj)
            doc.pop("_id")
            for f in doc.keys():
                f_str = "%s.**" % (f)
                idx_query = {f_str:1}
                idx_name  = {"name":"idx_on_%s" % (f)}
                dbh[coll].create_index([(f_str, pymongo.TEXT)], default_language='english')
                print idx_query, idx_name
                sys.exit()

    except pymongo.errors.ServerSelectionTimeoutError as err:
        print err
        return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
    except pymongo.errors.OperationFailure as err:
        print err
        return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

