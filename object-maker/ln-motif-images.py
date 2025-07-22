import os,sys
from optparse import OptionParser
import glob
import json
import subprocess
import datetime


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
    
    file_list = glob.glob(new_release_dir + "/jsondb/motifdb/*.json")
   
    for frmt in ["png", "svg", "json"] :
        d = new_release_dir + "/glycanimages_snfg_%s/" % (frmt)
        os.chdir(d)
        for in_file in file_list:
            doc = json.loads(open(in_file, "r").read())
            cmd = "ln -s %s.%s %s.%s" % (doc["glytoucan_ac"],frmt, doc["motif_ac"], frmt)
            x = subprocess.getoutput(cmd)

    os.chdir(".")





if __name__ == '__main__':
	main()

