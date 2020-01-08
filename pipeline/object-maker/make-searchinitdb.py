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
    
    res_obj = {
        "glycan_mass":{
            "native":{ "name":"Native", "min":0, "max":0 }
            ,"permethylated":{ "name":"Permethylated", "min":0, "max":0 }
        }
        ,"number_monosaccharides":{}, 
        "organism":[], 
        "glycan_type":[]
    }
    
    
    min_mass = 10000000.0
    max_mass = -min_mass
    min_mass_pme = 10000000.0
    max_mass_pme = -min_mass_pme

    min_monosaccharides = 100000
    max_monosaccharides = -min_monosaccharides

    seen = {"glycan_type":{}, "organism":{}}
    comp_dict = {}
    for obj in dbh[collection].find({}):
        glytoucan_ac = obj["glytoucan_ac"]
        if "mass" in obj:
            obj["mass"] = float(obj["mass"])
            if obj["mass"] > 0.0:
                min_mass = round(obj["mass"], 2) if obj["mass"] < min_mass else min_mass
                max_mass = round(obj["mass"], 2) if obj["mass"] > max_mass else max_mass
       
        if "mass_pme" in obj:
            obj["mass_pme"] = float(obj["mass_pme"])
            if obj["mass_pme"] > 0.0:
                min_mass_pme = round(obj["mass_pme"], 2) if obj["mass_pme"] < min_mass_pme else min_mass_pme
                max_mass_pme = round(obj["mass_pme"], 2) if obj["mass_pme"] > max_mass_pme else max_mass_pme


        if "number_monosaccharides" in obj:
            min_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] < min_monosaccharides else min_monosaccharides
            max_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] > max_monosaccharides else max_monosaccharides

        if "species" in obj:
            for o in obj["species"]:
                org_obj_str = json.dumps({ "name": o["name"], "id":o["taxid"]})
                seen["organism"][org_obj_str] = True
        if "composition" in obj:
            for o in obj["composition"]:
                res, res_name, res_count = o["residue"], o["name"], o["count"] 
                if res not in comp_dict:
                    comp_dict[res] = {"residue":res, "min":10000, "max":0 }
                if res_count < comp_dict[res]["min"]:
                    comp_dict[res]["min"] = res_count
                if res_count > comp_dict[res]["max"]:
                    comp_dict[res]["max"] = res_count

        if "classification" in obj:
            for o in obj["classification"]:
                type_name = o["type"]["name"]
                if type_name not in seen["glycan_type"]:
                    seen["glycan_type"][type_name] = {}
                if "subtype" in o:
                    subtype_name = o["subtype"]["name"]
                    if subtype_name not in seen["glycan_type"][type_name]:
                        seen["glycan_type"][type_name][subtype_name] = True 

    res_obj["glycan_mass"]["native"]["min"] = min_mass
    res_obj["glycan_mass"]["native"]["max"] = max_mass
    res_obj["glycan_mass"]["permethylated"]["min"] = min_mass_pme
    res_obj["glycan_mass"]["permethylated"]["max"] = max_mass_pme
    
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
    
    res_obj["composition"] = []
    for res in comp_dict:
        res_obj["composition"].append(comp_dict[res])

    return res_obj




def get_protein_search_init(dbh, collection):

    res_obj = {"protein_mass":{}, "organism":[]}
    min_mass = 10000000.0
    max_mass = -min_mass

    seen = {"protein_mass":{}, "organism":{}, "amino_acid":{}}
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
        if "glycosylation" in obj:
            for o in obj["glycosylation"]:
                seen["amino_acid"][o["residue"]] = True

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

    aa_dict = {
        "Arg":{ "name":"Arginine - Arg - R", "key":"R" },
        "Lys":{ "name":"Lysine - Lys - K", "key":"K" },
        "Asp":{ "name":"Aspartic acid - Asp - D", "key":"D" },
        "Glu":{ "name":"Glutamic acid - Glu - E", "key":"E" },
        "Gln":{ "name":"Glutamine - Gln - Q", "key":"Q" },
        "Asn":{ "name":"Asparagine - Asn - N", "key":"N" },
        "His":{ "name":"Histidine - His - H", "key":"H" },
        "Ser":{ "name":"Serine - Ser - S", "key":"S" },
        "Thr":{ "name":"Threonine - Thr - T", "key":"T" },
        "Tyr":{ "name":"Tyrosine - Tyr - Y", "key":"Y" },
        "Cys":{ "name":"Cysteine - Cys - C", "key":"C" },
        "Trp":{ "name":"Tryptophan - Trp - W", "key":"W" },
        "Ala":{ "name":"Alanine - Ala - A", "key":"A" },
        "Ile":{ "name":"Isoleucine - Ile - I", "key":"I" },
        "Leu":{ "name":"Leucine - Leu - L", "key":"L" },
        "Met":{ "name":"Methionine - Met - M", "key":"M" },
        "Phe":{ "name":"Phenylalanine - Phe - F", "key":"F" },
        "Val":{ "name":"Valine - Val - V", "key":"V" },
        "Pro":{ "name":"Proline - Pro - P", "key":"P" },
        "Gly":{ "name":"Glycine - Gly - G", "key":"G" }
    }
    res_obj["aa_list"] = []
    for aathree in seen["amino_acid"]:
        res_obj["aa_list"].append(aa_dict[aathree])

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
    
            out_file = path_obj["jsondbpath"] + "/searchinitdb/searchinitdb.json"
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

