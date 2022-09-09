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
        file_list = glob.glob("downloads/unicarbkb/*_glycosylation_current.csv")
        for in_file in file_list:
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                uniprotkb_ac = row[f_list.index("Protein")].split("/")[-1]
                glytoucan_ac = row[f_list.index("Toucan")] if "Toucan" in f_list else ""
                pos = row[f_list.index("Position")].split("^^")[0]
                pos = int(pos) if pos != "" else 0
                amino_acid = row[f_list.index("TypeAminoAcid")]
                amino_acid = amino_acid[0].upper() +  amino_acid[1:].lower()
                combo_id = "%s %s %s %s" % (uniprotkb_ac, pos, amino_acid, glytoucan_ac)
                seen[combo_id] = True
        nchecked, nfailed = 0, 0 
        for doc in dbh["c_protein"].find({}):
            uniprotkb_ac = doc["uniprot_canonical_ac"].split("-")[0]
            for obj in doc["glycosylation"]:
                glytoucan_ac = obj["glytoucan_ac"]
                if glytoucan_ac == "":
                    continue
                pos = str(obj["position"])
                amino_acid = obj["residue"]
                combo_id = "%s %s %s %s" % (uniprotkb_ac, pos, amino_acid, glytoucan_ac)
                nchecked += 1
                if combo_id not in seen:
                    nfailed += 1
                print "%s sites checked, %s sites failed" % (nchecked, nfailed) 

    except pymongo.errors.ServerSelectionTimeoutError as err:
        print err
        return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
    except pymongo.errors.OperationFailure as err:
        print err
        return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

