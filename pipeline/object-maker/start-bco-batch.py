import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
import datetime
import util
import hashlib

sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"





###############################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        data_dir = "reviewed/"

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

            

            seen = {"filename":{}, "bcoid":{}}
            for doc in dbh["c_bco"].find({}):
                doc.pop("_id")
                if "bco_id" in doc:
                    bco_id = doc["bco_id"]
                    bco_id_int = int(bco_id.split("_")[-1])
                    seen["bcoid"][bco_id] = True
                if "io_domain" in doc:
                    if "output_subdomain" in doc["io_domain"]:
                        for o in doc["io_domain"]["output_subdomain"]:
                            file_name = o["uri"]["filename"]
                            seen["filename"][file_name] = True

            name_dict = {
                "rykahsay@gwu.edu":"Robel Kahsay",
                "jeetvora@gwu.edu":"Jeet Vora",
                "rsn13@gwu.edu":"Rahi Navelkar",
                "xavierh@gwu.edu":"Xavier Holmes"
            }
            cont_list = []
            for e in name_dict:
                cont_list.append({
                    "name": name_dict[e],
                    "email":e,
                    "orcid": "",
                    "affiliation": "The George Washington University",
                    "contribution": ["createdBy"],
                })

            base_url = "http://data.glygen.org/datasets/reviewed/"
            in_file =  "tmp/ids.csv"
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            for row in data_frame["data"]:
                if row != []:
                    file_name = row[-1].strip()
                    desc = row[-2].strip()
                    if file_name in seen["filename"]:
                        print "BCO for %s exists" % (file_name)
                    else:
                        bco_obj = {
                            "usability_domain": [],
                            "provenance_domain": {
                                "name": desc,
                                "contributors": cont_list
                            },
                            "io_domain":{
                                "output_subdomain":[
                                    {
                                        "mediatype": "csv",
                                        "uri": {
                                            "filename": file_name,
                                            "access_time": "",
                                            "sha1_chksum": "",
                                            "uri": base_url + file_name
                                        }
                                    }
                                ]
                            },
                            "execution_domain": {
                                "software_prerequisites": [
                                    {
                                        "name": "Python", "version": "2.7.5", 
                                        "uri": {
                                            "uri": "https://www.python.org/download/releases/2.7.5/"
                                        }
                                    }
                                ],
                                "script":[
                                    {
                                        "uri":{
                                            "access_time":"2019-07-26T17:43:17-0500"
                                        }
                                    }
                                ]
                            }
                        }
                        seen["filename"][file_name] = True
                        bco_id = util.get_next_sequence_value(dbh["c_counters"], "bcoid")
                        bco_str_id = "000000"[0:6-len(str(bco_id))] + str(bco_id)
                        bco_url = config_obj["pathinfo"]["bcoprefix"] % (bco_str_id)
                        if bco_url in seen["bcoid"]:
                            print "BCO ID already exists %s " % (bco_url)
                            continue

                        bco_obj["bco_id"] = bco_url
                        bco_obj["provenance_domain"]["created"] = datetime.datetime.now().isoformat()
                        
                        bco_obj["checksum"] = hashlib.md5(json.dumps(bco_obj)).hexdigest()
                        print json.dumps(bco_obj, indent=4)
                        result = dbh["c_bco"].insert_one(bco_obj)
                        print "BCO created for %s" % (file_name)

        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

