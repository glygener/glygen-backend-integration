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
    (options,args) = parser.parse_args()
    for key in ([options.infile]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    in_file = options.infile
    in_json = json.loads(open(in_file, "r").read())

    outlines = ""
    for line in in_json["lines"]:
        outlines += "%s\n" % (line)

    with open("./requirements.txt", "w") as FW:
        FW.write("%s" % (outlines))


if __name__ == '__main__':
	main()
