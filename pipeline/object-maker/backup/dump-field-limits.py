import os,sys
import string
import commands
from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


def update_limit (obj, val):

    obj["min"] = val if val < obj["min"] else obj["min"]
    obj["max"] = val if val > obj["max"] else obj["max"]
    return




###############################
def main():


        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
	client = MongoClient('mongodb://localhost:27017')
	db = client[db_obj["dbname"]]


        limit_obj = {
            "glytoucan_ac":{"min":10000000000, "max":-1}
            ,"motif_name":{"min":10000000000, "max":-1}
            ,"uniprot_canonical_ac":{"min":10000000000, "max":-1}
            ,"uniprot_id":{"min":10000000000, "max":-1}
            ,"protein_name":{"min":10000000000, "max":-1}
            ,"gene_name":{"min":10000000000, "max":-1}
            ,"refseq_ac":{"min":10000000000, "max":-1}
            ,"pathway_id":{"min":10000000000, "max":-1}
            ,"pathway_name":{"min":10000000000, "max":-1}
            ,"disease_name":{"min":10000000000, "max":-1}
            ,"protein_sequence":{"min":10000000000, "max":-1}
        }
        

        for doc in db["c_glycan"].find({}):
            update_limit(limit_obj["glytoucan_ac"], len(doc["glytoucan_ac"]))
            for o in doc["motifs"]:
                update_limit(limit_obj["motif_name"], len(o["name"]))
        for doc in db["c_protein"].find({}):
            update_limit(limit_obj["uniprot_canonical_ac"], len(doc["uniprot_canonical_ac"]))
            update_limit(limit_obj["uniprot_id"], len(doc["uniprot_id"]))
            if "ac" in doc["refseq"]:
                update_limit(limit_obj["refseq_ac"], len(doc["refseq"]["ac"]))
            if "full" in doc["recommendedname"]:
                update_limit(limit_obj["protein_name"], len(doc["recommendedname"]["full"]))
            for o in doc["pathway"]:
                update_limit(limit_obj["pathway_name"], len(o["name"]))
                update_limit(limit_obj["pathway_id"], len(o["id"]))
            for o in doc["disease"]:
                update_limit(limit_obj["disease_name"], len(o["name"]))
            for o in doc["gene"]:
                update_limit(limit_obj["gene_name"], len(o["name"]))
            for o in doc["isoforms"]: 
                update_limit(limit_obj["protein_sequence"], len(o["sequence"]["sequence"]))

        limit_obj["enzyme_uniprot_canonical_ac"] = limit_obj["uniprot_canonical_ac"]

        print json.dumps(limit_obj, indent=4)


if __name__ == '__main__':
	main()

