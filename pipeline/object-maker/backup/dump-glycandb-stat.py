#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict


__version__="1.0"
__status__ = "Dev"


#######################################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

        jsondb_file = path_obj["jsondbpath"] + "/glycandb.json"
        out_obj_list = json.loads(open(jsondb_file, "r").read())
 	

	#key_list_1 = ['wurcs', 'glycoct', 'iupac']       
        #key_list_2 = ['mass']
        #key_list_3 = ['species', 'classification', 'crossref',  'glycosylation', 'enzyme']
        
        key_list_1 = []       
        key_list_2 = ['mass']
        key_list_3 = ['species', 'classification', 'glycosylation', 'enzyme']
        print ",".join(["glytoucan_id"] + key_list_1 + key_list_2 + key_list_3)
        for glycan_id in out_obj_list:
            obj = out_obj_list[glycan_id]
            if obj["species"] != []:
                print glycan_id

            #tv_list = [obj["id"]]
            #for key in key_list_1:
            #    tv_list.append("yes" if obj[key].strip() != "" else "no")
            #for key in key_list_2:
            #    tv_list.append("yes" if obj[key] > 0 else "no")
            #for key in key_list_3:
            #    tv_list.append("yes" if len(obj[key]) > 0 else "no")
            #print ",".join(tv_list)
            #print obj["id"], obj["glycosylation"]



if __name__ == '__main__':
        main()




