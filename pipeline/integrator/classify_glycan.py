import os,sys,glob
import json
import csv
from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON

__version__="1.1"
__status__ = "Dev"

###############################
class GLYCAN:
	def __init__(self, gtc_acc):
		self.gtc_acc = gtc_acc
		self.gtypeList = []
		self.gsubtypeList = []
		self.valid = True #will change to False if there is no type classification, then we don't print
	def add_type(self, gtype):
		self.gtypeList.append(gtype)
	def add_subtype(self, gsubtype):
		self.gsubtypeList.append(gsubtype)

	def determine_types(self):
		#type
		gtypeUnique = list(set(self.gtypeList)) #make unique
		if (len(gtypeUnique) == 0):
			self.valid = False #not N or O, don't print
		elif (len(gtypeUnique) == 1):
			self.gtype = gtypeUnique[0]
		else:
			print "ERROR: {} has non unique types {}".format(self.gtc_acc, ",".join(gtypeUnique))
			print "ERROR: {} has not been given a type and processing will fail"

		#subtype
		self.gsubtype = ""
		gsubtypeUnique = list(set(self.gsubtypeList))
		if (len(gsubtypeUnique) == 0):
			self.gsubtype = 'no subtype'
		elif (len(gsubtypeUnique) == 1):
			self.gsubtype = gsubtypeUnique[0]
		#this is a series of checks for known past issues
		#will return 'no subtype' for cases that cannot be solved by these 
		#better to provide no subtype than wrong subtype
		#for future analysis will print ids
		#for future additions beware of space before subtype name
		elif (len(gsubtypeUnique) == 2):
			#1 "N-Glycan high mannose","N-Glycan hybrid" are "N-Glycan hybrid" 
			if ((' high mannose' in gsubtypeUnique) and (' hybrid' in gsubtypeUnique)):
				self.gsubtype = ' hybrid'
		elif (len(gsubtypeUnique) == 3):
			#2 O-glycan core 1 + core 6 = core 2
			if ((" core 1" in gsubtypeUnique) and (" core 6" in gsubtypeUnique) and (" core 2" in gsubtypeUnique)):
				self.gsubtype = ' core 2'
			#2 O-glycan core 3 + core 6 = core 4
			if ((" core 3" in gsubtypeUnique) and (" core 6" in gsubtypeUnique) and (" core 4" in gsubtypeUnique)):
				self.gsubtype = ' core 4'
		if not(self.gsubtype):
			print "WARNING: {} has no discernable subtype. It contains motifs {}".format(self.gtc_acc, ",".join(gsubtypeUnique))
			self.gsubtype = 'no subtype'


	def typeRowList(self):
		return [self.gtc_acc, self.gtype]
	def subtypeRowList(self):
		return [self.gtc_acc, self.gsubtype]

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
	config_json = json.loads(open("/software/glygen/conf/config-1.json", "r").read())
	out_file1 = config_json["pathinfo"]["glytoucan_output_text"]+"/"+"glycan_type.txt"
	out_file2 = config_json["pathinfo"]["glytoucan_output_text"]+"/"+"glycan_subtype.txt"

	#NOTE these motif names must include TYPE (SPACE) Extra description or code will break
	motif_isSubtype = ["N-Glycan complex","N-Glycan high mannose","N-Glycan hybrid",
                           "O-Glycan core 1","O-Glycan core 2","O-Glycan core 3","O-Glycan core 4",
                           "O-Glycan core 5","O-Glycan core 6","O-Glycan core 7","O-Glycan core 8"]



	motif_notSubtype =["N-Glycan core basic","N-Glycan truncated motif. First GlcpNAC cut off",
                           "O-Glycan core 1 fuzzy","O-Glycan core 2 fuzzy","O-Glycan core 3 fuzzy",
                           "O-Glycan core 4 fuzzy","O-Glycan core 5 fuzzy","O-Glycan core 6 Fuzzy",
                           "O-Glycan core 7 fuzzy"]


	#TODO make GLYCAN object handle motifs directly
	glycan_dict = {}
	with open(in_file, "r") as MR:
		data = csv.reader(MR, delimiter='\t', quotechar='"')
		for row in data:
			line = row[0:3:2] #takes GlytoucanAccession Motif, one motif per line
			gtc_acc = line[0]
			motif = line[1]
			if not (gtc_acc in glycan_dict):
				glycan_dict[gtc_acc] = GLYCAN(gtc_acc)
			if motif in motif_isSubtype:
				glycan_dict[gtc_acc].add_type(motif.split(' ')[0])
				glycan_dict[gtc_acc].add_subtype(motif.split('-')[1].replace('Glycan',''))
			elif motif in motif_notSubtype:
				glycan_dict[gtc_acc].add_type(motif.split(' ')[0])
				#glycan_dict[gtc_acc].add_subtype('no subtype') #TODO is 'no subtype' necessary?, can we just not have an entry?

	

	#writing
	with open(out_file1, 'w') as csvfile:
		gtype_writer = csv.writer(csvfile, delimiter='\t')
		header_list = ["glytoucan_acc","type"]
		gtype_writer.writerow(header_list)
		with open(out_file2, 'w') as csvfile:
			gsubtype_writer = csv.writer(csvfile, delimiter='\t')
			header_list = ["glytoucan_acc","subtype"]
			gsubtype_writer.writerow(header_list)
			#Loop through dictionary, use class method to compute best guess subtype, write output
			#TODO will need to implement screening on glycan sequence info to handle edge cases, agreed
			#upon fuzzy fix
			for glycan in glycan_dict:
				glycan_dict[glycan].determine_types()
				if (glycan_dict[glycan].valid):
					gtype_writer.writerow(glycan_dict[glycan].typeRowList())
					gsubtype_writer.writerow(glycan_dict[glycan].subtypeRowList())

	
if __name__ == '__main__':
        main()
