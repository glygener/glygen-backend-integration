import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient



def load_species_dict(in_file):

    species2taxid = {}
    lines = open(in_file, "r").read().split("\n")[1:]
    for line in lines:
        if line.strip() == "":
            continue
        row = line[1:-1].split("\",\"")
        if row[-2] != "yes":
            continue
        tax_id = int(row[0])
        tax_name = row[2]
        if tax_name not in species2taxid:
            species2taxid[tax_name] = tax_id

    return species2taxid



def get_old_stat_obj(dbh, species_file):


    species2taxid = load_species_dict(species_file)
    stat_obj = {}
    for species in species2taxid:
        tax_id = species2taxid[species]
        tax_id_str = str(tax_id)
        stat_obj[tax_id_str] = {}
        stat_obj[tax_id_str]["species"] = species
        seen = {}
        for doc in dbh["c_glycan"].find({"species.taxid": {'$eq':tax_id}}):
            glytoucan_ac = doc["glytoucan_ac"]
            seen[glytoucan_ac] = True

        for doc in dbh["c_protein"].find({"$and":[{"species.taxid": {'$eq':tax_id}} , {"glycosylation": {'$gt':[]}}]}):
            for o in doc["glycosylation"]:
                glytoucan_ac = o["glytoucan_ac"]
                #if glytoucan_ac != "":
                #    seen[glytoucan_ac] = True
                if glytoucan_ac != "" and glytoucan_ac not in seen:
                    print glytoucan_ac, tax_id

        stat_obj[tax_id_str]["glycans"] = len(seen.keys())

    return stat_obj



def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
    db_obj = config_obj[config_obj["server"]]["dbinfo"]

    db_name = "glydb"
    ver = "1.9.6"

    misc_dir = "/data/shared/glygen/releases/data/v-%s/misc/" % (ver)
    user_info = {
        "glydb":{"mongodbuser":"glydbadmin" ,"mongodbpassword":"glydbpass"},
        "glydb_beta":{"mongodbuser":"glydb_betaadmin","mongodbpassword":"glydb_betapass"},
        "glydb_backup":{"mongodbuser":"glydb_backupadmin","mongodbpassword":"glydb_backuppass"}
    }

    client = pymongo.MongoClient('mongodb://localhost:27017',
        username=user_info[db_name]["mongodbuser"],
        password=user_info[db_name]["mongodbpassword"],
        authSource=db_name,
        authMechanism='SCRAM-SHA-1',
        serverSelectionTimeoutMS=10000
    )
    client.server_info()
    dbh = client[db_name]
        
    species_file = misc_dir + "/species_info.csv"
    old_stat_obj = get_old_stat_obj(dbh, species_file)
    print json.dumps(old_stat_obj, indent=4)


if __name__ == '__main__':
    main()

