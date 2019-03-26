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

        
        record_count = 0 
        out_obj_list = []
	for id1 in db_obj:
                

                out_obj_list.append(out_obj)
                record_count += 1
                if record_count%1000 == 0:
                    print " ... compiled %s objects" % (record_count)

        print " ... final compiled %s objects" % (record_count)

        record_count = 0
        fout_obj = {}
	for obj in out_obj_list:
	    cond_list= []
            cond_list.append(obj["mass"] > 0)
            cond_list.append(obj["number_monosaccharides"] > 0)
            if False not in cond_list:
                #remove empty value properties
                #clean_obj(obj)
                fout_obj[obj["glytoucan_ac"]] = order_obj(obj)
                record_count += 1
                if record_count%1000 == 0:
                    print " ... filtered %s objects" % (record_count)

        out_file = path_obj["jsondbpath"] + "glycandb.json"
        with open(out_file, "w") as FW:
	    FW.write("%s\n" % (json.dumps(fout_obj, indent=4)))


        print " ... final filtered in: %s objects" % (record_count)



if __name__ == '__main__':
        main()



