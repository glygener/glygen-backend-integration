import os,sys
from optparse import OptionParser
from rdflib import Graph
from lxml import etree
import json


__version__="1.0"
__status__ = "Dev"


###############################
def main():


	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--infile",action="store",dest="infile",help="Input file")
	parser.add_option("-d","--outdir",action="store",dest="outdir",help="Output directory")

	(options,args) = parser.parse_args()

	for file in ([options.infile, options.outdir]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	inFile = options.infile
	outDir = options.outdir

	FR = open(inFile, "r")
        i = 1
	
	seen1 = {}
	seen2 = {}
	badLines = []
	for line in FR:
                if "rdf:ID" in line:
                        val = line.split("rdf:ID=")[1].split("\"")[1]
			if val in seen1:
                                badLines.append(i)
			seen1[val] = 1
			if ":" in val:
                                seen2[val] = 1

		if i%1000000 == 0:
                        n = len(seen2.keys()) + len(badLines)
			print "parsed %s lines, found %s bad values" % (i, n)
                i += 1
	FR.close()


	fileName = os.path.basename(inFile)
	badValuesFile = outDir + "/" + fileName + "-badvalues"

	badValueList = seen2.keys()
	jsonObj = {"badlines":badLines, "badids":badValueList}
	with open(badValuesFile, "w") as FW:
		FW.write("%s\n" % (json.dumps(jsonObj, indent=4)))




if __name__ == '__main__':
	main()


