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

	fileName = os.path.basename(inFile)
	badValuesFile = outDir + "/" + fileName + "-badvalues"
	outFile = outDir + "/" + fileName

	jsonObj = json.loads(open(badValuesFile).read())
	FW = open(outFile, "w")	
	FR = open(inFile, "r")
	i = 1
	for line in FR:
		if "rdf:ID" in line:
			for w in jsonObj["badids"]:
				if line.find(w) != -1:
					ww = w.replace(":", "_")
					line = line.replace(w, ww)
                if i not in jsonObj["badlines"]:
			FW.write("%s" %(line))
		if i%1000000 == 0:
                        print "parsed %s lines" % (i)
		i += 1
	FR.close()
	FW.close()



if __name__ == '__main__':
	main()


