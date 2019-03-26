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


        tmpl_hash = {}
        with open("outdir/mapping.txt", 'r') as FR:
            data_frame = csv.reader(FR, delimiter='|', quotechar='"')
            rowCount = 0
            for row in data_frame:
                rowCount += 1
                tmpl_hash[row[0]] = row[1]


        jsondb_file = path_obj["jsondbpath"] + "/proteindb.json"
        out_obj_list = json.loads(open(jsondb_file, "r").read())
        
        for obj_id in out_obj_list:
            obj = out_obj_list[obj_id]
            for key1 in obj:
                if key1 == "mass":
                    triple = tmpl_hash[key1] % (obj_id, obj[key1])
                    print triple



if __name__ == '__main__':
        main()




