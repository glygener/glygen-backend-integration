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
    species_obj = config_obj["speciesinfo"]
    
    global path_obj
    path_obj = config_obj["pathinfo"]
    ds_obj = json.loads(open(path_obj["misc"] +  "dataset-masterlist.json", "r").read())
   

    primary_field_dict = {
            "protein":"uniprotkb_canonical_ac", 
            "glycan":"glytoucan_ac", 
            "proteoform":"uniprotkb_canonical_ac", 
            "miscprotein":"uniprotkb_canonical_ac", 
            "miscproteoform":"uniprotkb_canonical_ac"
    }


    obj_list = []
    for k_one in ["protein", "glycan", "proteoform", "miscprotein", "miscproteoform"]:
        obj = ds_obj[k_one]
        for k_two in obj:
            for ds in obj[k_two]:
                species_list = ["human", "mouse", "rat"]
                species_list = [k_two] if k_two != "common" else species_list
                primary_field = primary_field_dict[k_one]
                cat_obj = {"species":species_list, "molecule":k_one}
                o = {"name":ds, "format":"csv", "primaryfield":primary_field, "categories":cat_obj}
                obj_list.append(o)
    
    for ds in ds_obj["speciesagnostic"]:
        species_list = []
        primary_field = ""
        cat_obj = {"species":species_list, "molecule":""}
        o = {"name":ds, "format":"csv", "primaryfield":primary_field, "categories":cat_obj} 
        obj_list.append(o)
    
    
    print json.dumps(obj_list, indent=4)

    





if __name__ == '__main__':
	main()

