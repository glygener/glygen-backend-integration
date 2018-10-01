import os,sys
import json
import glob
import csv
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__="1.0"
__status__ = "Dev"

###############################
def main():

	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--infile",action="store",dest="infile",help="NT input file")
	parser.add_option("-d","--databasefile",action="store",dest="databasefile",help="cross_reference_list")
	(options,args) = parser.parse_args()
	for file in ([options.infile]):
		if not (file):
			parser.print_help()
			sys.exit(0)
	
	in_file = options.infile
	database_file =  options.databasefile
	config_json = json.loads(open("../conf/config-1.json", "r").read())

	if "human" in in_file:
		out_file = config_json["pathinfo"]["intermediate"]+"/"+"human_protein_database_urltemplate.csv"
	else:
		out_file = config_json["pathinfo"]["intermediate"]+"/"+"mouse_protein_database_urltemplate.csv"
    main_databases =[]
	with open(database_file, 'r') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
		rowcount = 0
		for row in csvreader:
			rowcount += 1
			if rowcount==1:
				continue
			else:
				main_databases += row

	url_database = {}
	databases = []
	data = []
	database_label =""
	with open(in_file, "r") as FR:
		for line in FR:
			row = line.split(' ')
			if "urlTemplate" in row[1]:
				database = row[0].split("/")[-1].replace(">","")
				database = database.lower()
				url_template = row[2].split('"')[1]
				if database not in url_database:
					url_database[database] = url_template
					

	header_list = ["database","url"]
	with open(out_file, 'w') as csvfile:  
		writer = csv.writer(csvfile)
		writer.writerow(header_list)
		for element in main_databases :
			row = [element]
			if element in url_database:	
				row.append(url_database[element])
			writer.writerow(row)

if __name__ == '__main__':
        main()
