import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient

sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"



###############################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        data_dir = "unreviewed/"

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


            field_dict = {}
            in_file = data_dir + "/field_names.csv"
            libgly.load_sheet_as_dict(field_dict, in_file, ",", "field_name_current")

            seen = {}
            for doc in dbh["c_bco"].find({}):
                bco_id = doc["bco_id"]
                file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"]
                file_ext = file_name.split(".")[-1]
                in_file = data_dir + "/" + file_name
                if os.path.isfile(in_file) == False:
                    continue
                
                if file_name not in seen:
                    seen[file_name] = True
                    if file_ext == "csv":
                        data_frame = {}
                        libgly.load_sheet(data_frame, in_file, ",")
                        stat_dict = {}
                        for row in data_frame["data"]:
                            for j in xrange(0, len(row)):
                                field = data_frame["fields"][j]
                                if field not in stat_dict:
                                    stat_dict[field] = []
                                stat_dict[field].append(row[j])
                        stat_obj = []
                        for field_name in stat_dict:
                            n = len(sorted(set(stat_dict[field_name])))
                            d = field_dict["data"][field_name][0][1] if field_name in field_dict["data"] else ""
                            stat_obj.append({"fieldname":field_name, "fielddesc":d, "uniquevalues":n})
                       

                        query_obj = {"bco_id":bco_id}
                        update_obj = {
                            "error_domain":{
                                "emperical_error": {
                                    "outputfields":stat_obj
                                },
                                "algorithmic_error":{}
                            }
                        }
                        result = dbh["c_bco"].update_one(query_obj, {'$set': update_obj}, upsert=True)
                        print "updated %s, %s" % (bco_id, file_name)

                

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

