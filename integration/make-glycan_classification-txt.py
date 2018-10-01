import os,sys,glob
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
    parser.add_option("-i","--infile",action="store",dest="infile",help="allmotif text file")
    (options,args) = parser.parse_args()
    for file in ([options.infile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

	in_file = options.infile
	config_json = json.loads(open("../conf/config-1.json", "r").read())
	out_file1 = config_json["pathinfo"]["intermediate"]+"/"+"glycan_type.text"
	out_file2 = config_json["pathinfo"]["intermediate"]+"/"+"glycan_subtype.txt"
	glycan_list = ["N-Glycan complex","N-Glycan high mannose","N-Glycan hybrid","O-Glycan core 1","O-Glycan core 2","O-Glycan core 3","O-Glycan core 4","O-Glycan core 5","O-Glycan core 6","O-Glycan core 7"]
	
	with open(out_file1, 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t')
		header_list = ["glytoucan_acc","type"]
		writer.writerow(header_list)
		with open(out_file2, 'w') as csvfile:
			writer1 = csv.writer(csvfile, delimiter='\t')
			header_list = ["glytoucan_acc","subtype"]
			writer1.writerow(header_list)
			with open(in_file, "r") as MR:
				data = csv.reader(MR, delimiter='\t', quotechar='"')
				for row in data:
					line = row[0:3:2]
					for item in glycan_list:
						if line[1] == item:
							glycan_type = [line[1].split(' ')[0]]
							glycan_subtype = [line[1].split('-')[1].replace('Glycan','')]
							gly_type = [line[0]]
							gly_type += glycan_type
							gly_subtype = [line[0]]
							gly_subtype += glycan_subtype
							writer.writerow(gly_type)
							writer1.writerow(gly_subtype)
	
if __name__ == '__main__':
        main()
