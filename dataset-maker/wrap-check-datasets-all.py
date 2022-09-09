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


    (options,args) = parser.parse_args()
    for file in ([options.molecule]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    molecule = options.molecule

    config_obj = json.loads(open("conf/config.json", "r").read())
    
    global path_obj
    path_obj = config_obj["pathinfo"]
    
    ds_obj_list = json.loads(open(path_obj["misc"] +  "dataset-masterlist.json", "r").read())     

    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if mol != molecule:
            continue
        if ds_name in ["homolog_alignments", "isoform_alignments"]:
            continue
        if ds_format in ["nt", "gif", "mp4"]:
            continue
        if obj["categories"]["species"] == []:
            out_file = path_obj["unreviewed"] + "%s_%s.%s" % (molecule, ds_name, ds_format)
            cmd = "wc %s" % (out_file)
            n_lines = commands.getoutput(cmd).strip().split(" ")[0]
            flag = "failed" if n_lines in ["0", "1"] or os.path.isfile(out_file) == False else "passed"
            print flag, "%s %s %s" % (molecule, ds_name, ds_format)
        else:
            for species in obj["categories"]["species"]:
                out_file = path_obj["unreviewed"] + "%s_%s_%s.%s" % (species,molecule,ds_name, ds_format)
                cmd = "wc %s" % (out_file)
                n_lines = commands.getoutput(cmd).strip().split(" ")[0]
                flag = "failed" if n_lines in ["0", "1"] or os.path.isfile(out_file) == False else "passed"
                print flag, "%s %s %s %s" % (species,molecule,ds_name, ds_format)



if __name__ == '__main__':
	main()

