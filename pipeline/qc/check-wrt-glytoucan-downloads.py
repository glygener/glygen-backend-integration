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




###############################
def main():

    config_obj = json.loads(open("conf/config.json", "r").read())
    db_obj = config_obj[config_obj["server"]]["dbinfo"]

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
        
        seen = {}
        file_list = glob.glob("downloads/glytoucan/current/export/allglycan.tsv")
        for in_file in file_list:
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
                seen[glytoucan_ac] = True

        nchecked, nfailed = 0, 0 
        for doc in dbh["c_glycan"].find({}):
            glytoucan_ac = doc["glytoucan_ac"]
            nchecked += 1
            if glytoucan_ac not in seen:
                nfailed += 1
            print "%s checked, %s failed" % (nchecked, nfailed) 

    except pymongo.errors.ServerSelectionTimeoutError as err:
        print err
        return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
    except pymongo.errors.OperationFailure as err:
        print err
        return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

