import os,sys
import string
import json
from optparse import OptionParser
import csv

__version__="1.0"
__status__ = "Dev"


###############################
def main():

	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--infile",action="store",dest="infile",help="CSV input file")
	parser.add_option("-x","--excludelist",action="store",dest="excludelist",help="Field index list to be excluded '[2,3]'")
	parser.add_option("-o","--organism",action="store",dest="organism",help="human or mouse")
	(options,args) = parser.parse_args()
	for file in ([options.infile, options.excludelist, options.organism]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	in_file = options.infile
	exclude_list = json.loads(options.excludelist)
	in_filename = os.path.basename(in_file)
	in_filename, in_fileext = os.path.splitext(os.path.basename(in_file))
	organism = options.organism
	config_json = json.loads(open("../conf/config-1.json", "r").read())
	out_file = config_json["pathinfo"]["intermediate"] +"/" + organism +"_"+ in_filename + ".csv"
	log_file = config_json["pathinfo"]["intermediate"] +"/" + organism +"_" + in_filename + ".log"
	
	ac2canonical = {}
	idmapfile = config_json["pathinfo"]["reviewed"] +"/human_protein_idmapping.csv"
	with open(idmapfile, 'r') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
		rowCount = 0
		for row in csvreader:
			rowCount += 1
			if rowCount == 1:
				continue
			ac = row[0].split("-")[0]
			ac2canonical[ac] = row[0]
	fw2 = open(log_file, "w")
	with open(outFile, 'w') as output_csv:
		csvwriter = csv.writer(output_csv)
		count = {}
		field_list = []
		with open(in_file, 'r') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
			rowCount = 0
			for row in csvreader:
				rowCount += 1
				if rowCount == 1:
					line = [row[0]]				
					for i in xrange(1, len(row)):
						if i not in exclude_list:
							line += [row[i]]
					csvwriter.writerow(line)
				else:
					ac = row[0].strip()
					if ac not in ac2canonical:
						fw2.write("%s,entry does NOT exist anymore\n" % (ac))
					else:
						line = [ac2canonical[ac]]
						for i in xrange(1, len(row)):
							if i not in exclude_list:
								line += [row[i]]
						csvwriter.writerow(line)							
	fw2.close()

	print "\n\tOutput file: %s" % (out_file)
	print "\tLog file: %s\n\n" % (log_file)
	
if __name__ == '__main__':
        main()

