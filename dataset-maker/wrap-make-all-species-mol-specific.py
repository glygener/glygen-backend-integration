import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient

import libgly



__version__="1.0"
__status__ = "Dev"



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-m","--molecule",action="store",dest="molecule",help="glycan/protein/proteoform")
    parser.add_option("-s","--species",action="store",dest="species",help="species (optional)")


    (options,args) = parser.parse_args()
    for file in ([options.molecule]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    molecule = options.molecule
    species_name = options.species
    
    config_obj = json.loads(open("conf/config.json", "r").read())
    
    global path_obj
    path_obj = config_obj["pathinfo"]
    ds_obj_list = json.loads(open(path_obj["misc"] +  "dataset-masterlist.json", "r").read())     


    cmd_list_one, cmd_list_two = [], []
    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if mol != molecule:
            continue
        if obj["categories"]["species"] != []:
            for species in obj["categories"]["species"]:
                if species_name != None and species != species_name:
                    continue
                out_file = path_obj["unreviewed"] + "%s_%s_%s.%s" % (species, molecule, ds_name, ds_format)
                cmd = "/usr/bin/python make-%s-dataset.py -s %s -d %s " % (molecule,species, ds_name)
                if ds_name not in ["isoform_alignments", "ntdata"]:
                    cmd += " > %s" % (out_file)
                if ds_name.find("citations_") != -1:
                    cmd_list_two.append(cmd)
                else:
                    cmd_list_one.append(cmd)

    for cmd in cmd_list_one + cmd_list_two:
        x = commands.getoutput(cmd)
        #print cmd






if __name__ == '__main__':
	main()

