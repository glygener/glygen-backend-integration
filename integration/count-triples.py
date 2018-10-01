import os,sys
import json
import csv


from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


__version__="1.0"
__status__ = "Dev"




###############################
def main():

	usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i","--infile",action="store",dest="infile",help="Input n3 file")
        parser.add_option("-d","--domain",action="store",dest="domain",help="Domain (up:Protein)")
        parser.add_option("-p","--predicate",action="store",dest="predicate",help="Predicate (up:sequence)")


        (options,args) = parser.parse_args()
        for file in ([options.infile, options.domain, options.predicate]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

        prefixmap = json.loads(open("./prefixmap.json", "r").read())
	infile = options.infile
        domain = options.domain
        predicate = options.predicate

        domain_prefix = prefixmap["domain"][domain]
        predicate_prefix = prefixmap["predicate"][predicate]


        n = 0
        with open(infile, "r") as FR:
            for line in FR:
                if len(line.strip()) == 0:
                    continue
                parts = line.strip().split(" ")
                c1 = parts[0].find(domain_prefix) != -1
                c2 = parts[0].find("#") == -1
                c3 = parts[1].find(predicate_prefix) != -1
                if c1 and c2 and c3:
                    n += 1

        print "%s --> %s (%s triples)"  % (domain, predicate, n)


if __name__ == '__main__':
	main()
