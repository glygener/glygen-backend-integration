import os,sys
import string
import subprocess
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
import subprocess
import datetime
#import requests
from Bio import SeqIO


__version__="1.0"
__status__ = "Dev"



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-v","--ver",action="store",dest="ver",help="Version")

    
    (options,args) = parser.parse_args()
    for key in ([options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    ver  = options.ver
    data_release_dir = "/data/shared/glygen/releases/data/"
    new_release_dir = data_release_dir + "v-%s/" % (ver)
    os.chdir(new_release_dir + "/glycanimages_snfg/")
    
    file_list = glob.glob(new_release_dir + "/jsondb/motifdb/*.json")
    for in_file in file_list:
        doc = json.loads(open(in_file, "r").read())
        cmd = "ln -s %s.png %s.png" % (doc["glytoucan_ac"],doc["motif_ac"])
        x = subprocess.getoutput(cmd)
        #print (cmd)

    os.chdir(".")





if __name__ == '__main__':
	main()

