import os,sys,glob
import json
import csv
from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON

__version__="1.1"
__status__ = "Dev"

"""
BCF Nov. 21, 2018
v 1.0 took first human orthologue accession number from infile1 and last human gene symbol from infile1
this version fixes that and ensures that the accession number and gene symbol match in infile1
TODO add preference for human gene symbol that matches mouse gene symbol and add commenting
"""
###############################
def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--infile1",action="store",dest="infile1",help="mouse protein ortholog csv file")
	parser.add_option("-d","--infile2",action="store",dest="infile2",help="gene info for humans")
	parser.add_option("-s","--infile3",action="store",dest="infile3",help="mouse glycan enzyme ortholog csv file")
	(options,args) = parser.parse_args()
	for file in ([options.infile1, options.infile2, options.infile3]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	in_file1 = options.infile1
	in_file2 = options.infile2
	in_file3 = options.infile3

	config_json = json.loads(open("/software/glygen/conf/config-1.json", "r").read())
	out_file = config_json["pathinfo"]["glytoucan_output_csv"]+"/"+"human_glycan_enzyme.csv"
	
	accession = {}
	gene_sym = {}
	with open(in_file1, "r") as FR:
		dataGrid = csv.reader(FR, delimiter=',', quotechar='"')
		rowCount = 0
		for row in dataGrid:
			cannonical = row[0]
			#cannonical = row[0].split('-')[0]
			#print cannonical
			ortholog_cannonical = row[3]
			gene = row[5]
			gene_symbol = row[4]
			if cannonical not in accession:
				accession[cannonical] = []
				gene_sym[cannonical] = []
			accession[cannonical].append(ortholog_cannonical)
			accession[cannonical].append(gene)
			gene_sym[cannonical].append(gene_symbol)
	gene_id = {}
	with open(in_file2, "r") as GR:
		data = csv.reader(GR, delimiter='\t', quotechar='"')
		rowCount = 0
		for row in data:
			gene_id[row[1]] = row[8]
	cann_acc = []
	header_list = ["glytoucan_acc","uniprotkb_acc_canonical_enzyme","gene_symbol_enzyme","gene_id_enzyme","enzyme_name","tax_id_enzyme","organism_name_enzyme"]
	with open(out_file, 'w') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(header_list)
		with open (in_file3,"r") as CR:
			datagrid = csv.reader(CR, delimiter=',', quotechar='"')
			rowCount = 0
			for row in datagrid:
				acc = row[1]
				line = [row[0]]
				if acc in accession:
					line.append(accession[acc][0])
					line.append(gene_sym[acc][0])
					line.append(accession[acc][1])
					if accession[acc][1] in gene_id:
						line += [gene_id[accession[acc][1]]]	
						line += ["9606"]
						line += ["Homo sapiens"]
						#print line
						writer.writerow(line)	

if __name__ == '__main__':
        main()
