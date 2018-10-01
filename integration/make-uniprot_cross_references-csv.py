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
	parser.add_option("-i","--infile",action="store",dest="infile",help="cross reference n3 file")
	parser.add_option("-d","--databasefile",action="store",dest="databasefile",help="cross reference list")
	(options,args) = parser.parse_args()
	for file in ([options.infile, options.databasefile]):
		if not (file):
			parser.print_help()
			sys.exit(0)
    
	in_file = options.infile
	database_file = options.databasefile
	config_json = json.loads(open("../conf/config-1.json", "r").read())
	if "human" in in_file:
		out_file = config_json["pathinfo"]["intermediate"]+"/"+"human_protein_uniprot_crossreference.csv"
	else:
		out_file = config_json["pathinfo"]["intermediate"]+"/"+"mouse_protein_uniprot_crossreference.csv"

	main_databases =[]
	with open(database_file, 'r') as csvfile:
		csvReader = csv.reader(csvfile, delimiter=',', quotechar='|')
		rowcount = 0
		for row in csvReader:
			rowcount += 1
			if rowcount==1:
				continue
			else:
				main_databases += row
	
	uniprotacc_database = {}
	comment_database = {}
	databases = []
	data = []
	database_label =""
	rowcount = 0
	with open(in_file, "r") as FR:
		for line in FR:
			row = line.split(' ')
			object_part = row[2].split("/")
			if "comment" in line:
				row = line.split('"')
				database_label = row[1]
				database = row[0].split("/")[3]	
				d_id = row[0].split("/")[4].replace(">","").split(" ")[0]
				if database not in data:
					data += [database]
				if database not in comment_database:
					comment_database[database] = {}
				if d_id not in comment_database[database]:
					comment_database[database][d_id] = []
				if database_label not in comment_database[database][d_id]:
					comment_database[database][d_id].append(database_label)
			for item in object_part:
				if item in main_databases or "nextprot" in item:
					database = item.split(".")[0]
					uniprot_acc = row[0].split("/")[-1].replace(">","").split("#")[0]
					database_id = row[2].split("/")[-1].replace(">","")
					if database not in databases:
						databases += [database]
					if uniprot_acc not in uniprotacc_database:
						uniprotacc_database[uniprot_acc]={}
					if database not in uniprotacc_database[uniprot_acc]:
						uniprotacc_database[uniprot_acc][database]= []
					uniprotacc_database[uniprot_acc][database].append(database_id)
					
	rowcount = 0
	header_list = []
	for element in databases:
		rowcount += 1
		header_list += [element+"_id",element+"_label"]
	header_list = ["uniprotkb_acc"] + header_list
	with open(out_file, 'w') as csvfile:  
		writer = csv.writer(csvfile)
		writer.writerow(header_list)
		for uniprot_acc in uniprotacc_database:
			row = [uniprot_acc]
			for element in databases:
				if element in uniprotacc_database[uniprot_acc]:
					database_id = uniprotacc_database[uniprot_acc][element]	
					row.append("|".join(uniprotacc_database[uniprot_acc][element]))							
					if element in comment_database:
						temp=[]
						for item in database_id:
							if item in comment_database[element]:
								temp += comment_database[element][item]
						row.append("|".join(temp))
					else:
						row.append("")	
				else:
					row.append("")
					row.append("")
			writer.writerow(row)

if __name__ == '__main__':
        main()
