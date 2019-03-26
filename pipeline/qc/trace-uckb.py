import os,sys
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

from Bio import SeqIO
import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"


###############################
def main():

		
	usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

        (options,args) = parser.parse_args()
        for file in ([options.species]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

        
        species = options.species
        tax_id = "9606" if species == "human" else "10090"
       
  
        in_file = "/data/projects/glygen/downloads/uckb/human_glycosylation.csv"
        if tax_id == "10090":
            in_file = "/data/projects/glygen/downloads/uckb/mouse_glycosylation.csv"

        t_list = {}
	src_grid = []
        with open(in_file, 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                if row[-2][0:5] == "comp_":
                    continue
                row_count += 1
                if row_count == 1:
                    field_list = row
                else:
                    row_obj = {}
                    for j in xrange(0,len(field_list)):
                        field_name  = field_list[j]
                        row_obj[field_name] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                    src_grid.append(row_obj)
                    combo_id = row_obj["Protein"][0].split("/")[-1] + "," + row_obj["Position"][0].split("^^")[0]
                    if combo_id not in t_list:
                        t_list[combo_id] = []
                    t_list[combo_id] += row_obj["Id"]

        #for combo_id in t_list:
        #    print len(sorted(set(t_list[combo_id]))), combo_id, t_list[combo_id]
        #sys.exit()

        seen = {}
        for obj in src_grid:
            ac = obj["Protein"][0].split("/")[-1]
            pos = obj["Position"][0].split("^^")[0]
            amino_acid = obj["TypeAminoAcid"][0]
            uckb_id = obj["Id"][0]
            combo_id = "%s,%s,%s,%s" % (ac,pos,amino_acid,uckb_id) 
            seen[combo_id] = True



	client = MongoClient('mongodb://localhost:27017')
	db = client["glygen"]
        
        cond_list = []
        cond_list.append({"species.taxid": {'$eq': int(tax_id)}})
        cond_list.append({"glycosylation.evidence.database": {'$eq':"UniCarbKB"}})
        query_obj = {} if cond_list == [] else { "$and": cond_list }
        for doc in db["c_protein"].find(query_obj):
            canon = doc["uniprot_canonical_ac"]
            ac = canon.split("-")[0]
            for obj in doc["glycosylation"]:
                pos,amino_acid = obj["position"],obj["residue"]
                for o in obj["evidence"]:
                    if o["database"] == "UniCarbKB":
                        uckb_id = o["id"]
                        combo_id = "%s,%s,%s,%s" % (ac,pos,amino_acid,uckb_id)
                        print combo_id, combo_id in seen



if __name__ == '__main__':
	main()


