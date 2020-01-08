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


        dir_name = "glygen-beta" if db_name == "glydb_beta" else "glygen" 
        db_path = "/data/projects/%s/wwwdata/jsondb/" % (dir_name) 
        db_path = "jsondb/"

        
        coll_info = {
            "c_glycan":{
                "dbfile":db_path + "glycandb/*.json"
                ,"primaryid":"glycanid"
                ,"dbtype":1
            }
            ,"c_protein":{
                "dbfile":db_path + "proteindb/*.json"
                ,"primaryid":"proteinid"
                ,"dbtype":1
            }
            ,"c_alignment":{
                "dbfile":db_path + "alignmentdb/*.json"
                ,"primaryid":"alignmentid"
                ,"dbtype":1
            }
            ,"c_cluster":{
                "dbfile":db_path + "clusterdb/*.json"
                ,"primaryid":"clusterid"
                ,"dbtype":1
            }
            ,"c_searchinit":{
                "dbfile":db_path + "searchinitdb/searchinitdb.json"
                ,"primaryid":"searchinitid"
                ,"dbtype":2
            }
        }
        json_db_file = coll_info[coll]["dbfile"]
        primary_id = coll_info[coll]["primaryid"]

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

            flag = True
            reset_sequence_value(dbh, primary_id)
            nrecords = 0
            if coll_info[coll]["dbtype"] == 1:
                file_list = glob.glob(coll_info[coll]["dbfile"])
                for in_file in file_list:
                    doc = json.loads(open(in_file, "r").read())
                    doc["_id"]  = get_next_sequence_value(dbh, primary_id)
                    result = dbh[coll].insert_one(doc)	
                    nrecords += 1
            elif coll_info[coll]["dbtype"] == 2:
                in_file = coll_info[coll]["dbfile"]
                doc = json.loads(open(in_file, "r").read())
                result = dbh[coll].insert_one(doc)
                nrecords += 1
            print "Done! Loaded %s objects to %s " % (nrecords, coll)
        except pymongo.errors.ServerSelectionTimeoutError as err:
            print err
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            print err
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

