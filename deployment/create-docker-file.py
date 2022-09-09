import os,sys
import string
from optparse import OptionParser
import glob
import json


__version__="1.0"
__status__ = "Dev"



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-i","--infile",action="store",dest="infile",help="Input File")
    parser.add_option("-m","--mod",action="store",dest="mod",help="Module [glygen/argosdb]")
    parser.add_option("-s","--ser",action="store",dest="ser",help="Server [dev/tst/beta/prd]")
    parser.add_option("-v","--ver",action="store",dest="ver",help="Version")

    (options,args) = parser.parse_args()
    for key in ([options.infile, options.mod, options.ser, options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    in_file = options.infile
    mod = options.mod
    ser = options.ser
    ver = options.ver
    in_json = json.loads(open(in_file, "r").read())

    param_dict = {}
    for param in in_json["parameters"]:
        param_dict[param] = in_json["parameters"][param]
    
    param_dict["$mod"] = mod
    param_dict["$ser"] = ser
    param_dict["$ver"] = ver

    outlines = ""
    for line in in_json["lines"]:
        for param in param_dict:
            line = line.replace(param, param_dict[param])
        outlines += "%s\n\n" % (line)

    with open("./Dockerfile", "w") as FW:
        FW.write("%s" % (outlines))


if __name__ == '__main__':
	main()
