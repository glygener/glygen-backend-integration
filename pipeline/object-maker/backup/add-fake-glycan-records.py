#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict
from random import randint

__version__="1.0"
__status__ = "Dev"



#####################################
def order_obj(jsonObj):

    ordrHash = {"glytoucan_ac":1, "mass":2, "iupac":3, "wurcs":4, "glycoct":5,
            "species":6, "classification":7,"glycosylation":8,
            "enzyme":9, "crossref":10}

    for k1 in jsonObj:
        ordrHash[k1] = ordrHash[k1] if k1 in ordrHash else 1000
        if type(jsonObj[k1]) is dict:
            for k2 in jsonObj[k1]:
                ordrHash[k2] = ordrHash[k2] if k2 in ordrHash else 1000
                if type(jsonObj[k1][k2]) is dict:
                    for k3 in jsonObj[k1][k2]:
                        ordrHash[k3] = ordrHash[k3] if k3 in ordrHash else 1000
                    jsonObj[k1][k2] = OrderedDict(sorted(jsonObj[k1][k2].items(),
                        key=lambda x: float(ordrHash.get(x[0]))))
                elif type(jsonObj[k1][k2]) is list:
                    for j in xrange(0, len(jsonObj[k1][k2])):
                        if type(jsonObj[k1][k2][j]) is dict:
                            for k3 in jsonObj[k1][k2][j]:
                                ordrHash[k3] = ordrHash[k3] if k3 in ordrHash else 1000
                                jsonObj[k1][k2][j] = OrderedDict(sorted(jsonObj[k1][k2][j].items(), 
                                    key=lambda x: float(ordrHash.get(x[0]))))
            jsonObj[k1] = OrderedDict(sorted(jsonObj[k1].items(),
                key=lambda x: float(ordrHash.get(x[0]))))

    return OrderedDict(sorted(jsonObj.items(), key=lambda x: float(ordrHash.get(x[0]))))



#######################################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

        in_file = path_obj["jsondbpath"] + "glycandb.json"
                
        glycandb_objs = json.loads(open(in_file, "r").read())
        example_obj = json.loads(open("tmp/detail_example.json", "r").read())
            
        key_list_1 = ['wurcs', 'glycoct', 'iupac']       
        key_list_2 = ['mass', 'number_monosaccharides']
        key_list_3 = ['motifs', 'species',  'crossref',  'glycoprotein', 'enzyme']
        
        pool = {}
        for key in key_list_1 + key_list_2 + key_list_3:
            pool[key] = []

        for glycan_id in glycandb_objs:
            obj = glycandb_objs[glycan_id]
            for key in key_list_1:
                if obj[key].strip() != "":
                    pool[key].append(obj[key])
            for key in key_list_2:
                if obj[key] > 0:
                    pool[key].append(obj[key])
            for key in key_list_3:
                if len(obj[key]) > 0:
                    pool[key].append(obj[key])

        new_obj_list = {}
        for glycan_id in glycandb_objs:
            obj = glycandb_objs[glycan_id]
            if len(obj["glycoprotein"]) > 0:
                new_id = "GLY" + glycan_id[3:]
                new_obj_list[new_id] = {"glytoucan_ac":new_id}
                for key in key_list_1 + key_list_2 + key_list_3:
                    r = randint(0, len(pool[key])-1)
                    new_obj_list[new_id][key] = pool[key][r]
                new_obj_list[new_id]["classification"] = example_obj["classification"] 
        
        for new_id in new_obj_list:
            glycandb_objs[new_id] = order_obj(new_obj_list[new_id])


        for i in xrange(1, 10):
            new_id = "GLYMIN0" + str(i)
            r = randint(0, len(pool["wurcs"])-1)
            wurcs_id = pool["wurcs"][r]
            glycandb_objs[new_id] = {"glytoucan_ac":new_id, "wurcs":wurcs_id}

       
        out_file = path_obj["jsondbpath"] + "glycandbplus.json"
        with open(out_file, "w") as FW: 
            FW.write("%s\n" % (json.dumps(glycandb_objs, indent=4, sort_keys=True)))



if __name__ == '__main__':
        main()




