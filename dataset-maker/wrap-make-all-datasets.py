import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json

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

    filtered_obj_list_zero = []
    filtered_obj_list_one, filtered_obj_list_two, filtered_obj_list_three = [], [], []
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
                o = {"name":ds_name, "format":ds_format, "species":species}
                if ds_name.find("masterlist") != -1:
                    filtered_obj_list_zero.append(o)
                elif ds_name.find("canonicalsequences") != -1 or ds_name.find("allsequences") != -1:
                    filtered_obj_list_one.append(o)
                elif ds_name.find("citations_") != -1:
                    filtered_obj_list_three.append(o)
                else:
                    filtered_obj_list_two.append(o)
        else:
            o = {"name":ds_name, "format":ds_format}
            if molecule == "glycan" and ds_name == "masterlist":
                filtered_obj_list_one.append(o)
            elif ds_name.find("citations_") == -1:
                filtered_obj_list_two.append(o)
            else:
                filtered_obj_list_three.append(o)

 
    for o in filtered_obj_list_zero + filtered_obj_list_one + filtered_obj_list_two + filtered_obj_list_three:
        ds_name = o["name"]
        ds_format = o["format"]
        species = o["species"] if "species" in o else ""
        out_file = "unreviewed/%s_%s.%s" % (molecule, ds_name, ds_format)
        cmd = "/usr/bin/python make-species-agnostic-dataset.py -d %s " % (ds_name)
        if species != "":
            out_file = "unreviewed/%s_%s_%s.%s" % (species, molecule, ds_name, ds_format)
            cmd = "/usr/bin/python make-%s-dataset.py -s %s -d %s " % (molecule,species, ds_name)
        if ds_name not in ["isoform_alignments", "ntdata", "homolog_alignments", "glycan_images"]:
            cmd += " > %s" % (out_file)
        x = commands.getoutput(cmd)
        #print cmd





if __name__ == '__main__':
	main()

