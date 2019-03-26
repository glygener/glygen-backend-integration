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


def get_glycan_search_init(dbh, collection):
    
    res_obj = {"glycan_mass":{}, "number_monosaccharides":{}, 
                "organism":[], "glycan_type":[]}
    min_mass = 10000000.0
    max_mass = -min_mass
    min_monosaccharides = 100000
    max_monosaccharides = -min_monosaccharides

    seen = {"glycan_type":{}, "organism":{}}
    for obj in dbh[collection].find({}):
        glytoucan_ac = obj["glytoucan_ac"]
        if "mass" in obj:
            obj["mass"] = float(obj["mass"])
            min_mass = round(obj["mass"], 2) if obj["mass"] < min_mass else min_mass
            max_mass = round(obj["mass"], 2) if obj["mass"] > max_mass else max_mass
        
        if "number_monosaccharides" in obj:
            min_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] < min_monosaccharides else min_monosaccharides
            max_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] > max_monosaccharides else max_monosaccharides

        if "species" in obj:
            for o in obj["species"]:
                org_obj_str = json.dumps({ "name": o["name"], "id":o["taxid"]})
                seen["organism"][org_obj_str] = True
        
        if "classification" in obj:
            for o in obj["classification"]:
                type_name = o["type"]["name"]
                if type_name not in seen["glycan_type"]:
                    seen["glycan_type"][type_name] = {}
                if "subtype" in o:
                    subtype_name = o["subtype"]["name"]
                    if subtype_name not in seen["glycan_type"][type_name]:
                        seen["glycan_type"][type_name][subtype_name] = True 

    res_obj["glycan_mass"]["min"] = min_mass
    res_obj["glycan_mass"]["max"] = max_mass
    res_obj["number_monosaccharides"]["min"] = min_monosaccharides
    res_obj["number_monosaccharides"]["max"] = max_monosaccharides 

    for type_name in seen["glycan_type"]:
        o = {"name":type_name,"subtype":[]}
        for subtype_name in seen["glycan_type"][type_name]:
            o["subtype"].append(subtype_name)
        res_obj["glycan_type"].append(o)
    
    for org_obj_str in seen["organism"]:
        res_obj["organism"].append(json.loads(org_obj_str))


    res_obj["simple_search_category"] = [
         {"id":"any", "display":"Any"}
        ,{"id":"glycan", "display":"Glycan"}
        ,{"id":"protein", "display":"Protein"}
        ,{"id":"enzyme", "display":"Enzyme"}
        ,{"id":"organism", "display":"Organism"}
    ]
    return res_obj




def get_protein_search_init(dbh, collection):

    res_obj = {"protein_mass":{}, "organism":[]}
    min_mass = 10000000.0
    max_mass = -min_mass

    seen = {"protein_mass":{}, "organism":{}}
    for obj in dbh[collection].find({}):
        uniprot_canonical_ac = obj["uniprot_canonical_ac"]
        if "mass" in obj:
            mass = float(obj["mass"]["chemical_mass"])
            if mass == -1.0:
                continue
            min_mass = round(mass, 2) if mass < min_mass else min_mass
            max_mass = round(mass, 2) if mass > max_mass else max_mass
        if "species" in obj:
            for o in obj["species"]:
                org_obj_str = json.dumps({ "name": o["name"], "id":o["taxid"]})
                seen["organism"][org_obj_str] = True

    res_obj["protein_mass"]["min"] = min_mass
    res_obj["protein_mass"]["max"] = max_mass
    for org_obj_str in seen["organism"]:
        res_obj["organism"].append(json.loads(org_obj_str))

    res_obj["simple_search_category"] = [
         {"id":"any", "display":"Any"}
        ,{"id":"glycan", "display":"Glycan"}
        ,{"id":"protein", "display":"Protein"}
        ,{"id":"organism", "display":"Organism"}
        ,{"id":"disease", "display":"Disease"}
        ,{"id":"pathway", "display":"Pathway"}
    ] 

    return res_obj




###############################
def main():

        config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
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


            out_obj = {}
            out_obj["protein"] = get_protein_search_init(dbh, "c_protein")
            out_obj["glycan"] = get_glycan_search_init(dbh, "c_glycan")
    
            out_file = path_obj["jsondbpath"] + "searchinitdb.json"
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(out_obj, indent=4)))
        except pymongo.errors.ServerSelectionTimeoutError as err:
            print err
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            print err
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}




if __name__ == '__main__':
	main()

