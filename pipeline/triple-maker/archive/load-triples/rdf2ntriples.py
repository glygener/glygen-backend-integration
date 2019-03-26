import os,sys
from optparse import OptionParser
#from rdflib import Graph
#from lxml import etree
import json
import commands


__version__="1.0"
__status__ = "Dev"



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--infile",action="store",dest="infile",help="Input file")
    parser.add_option("-o","--outfile",action="store",dest="outfile",help="Output file")

    (options,args) = parser.parse_args()
    for file in ([options.infile, options.outfile]):
	if not (file):
	    parser.print_help()
	    sys.exit(0)
    
    in_file = options.infile
    out_file = options.outfile

    in_filename = os.path.basename(in_file)
    out_filename =  ('.').join(in_filename.split('.')[:-1])
	
    log_file = out_file + ".log"
    cmd = "/usr/local/bin/rapper -o ntriples " + in_file + " 1> " + out_file + " 2> " + log_file
    x = commands.getoutput(cmd)


if __name__ == '__main__':
    main()


