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
       
   

        in_file = "/data/projects/glygen/downloads/uniprot/uniprot-proteome-9606-homo-sapiens-2017-11.nt"
        if tax_id == "10090":
            in_file = "/data/projects/glygen/downloads/uniprot/uniprot-proteome-10090-mus-musculus-2017-11.nt"

        is_canon = {}
        ac2isoformlist = {}
        with open(in_file, "r") as FR:
            for line in FR:
                prd = "<http://purl.uniprot.org/core/sequence>"
                if line.find(prd) != -1:
                    ac = line.strip().split(" ")[0].split("/")[-1][:-1]
                    isoform = line.strip().split(" ")[-2].split("/")[-1][:-1]
                    if ac.find("#") != -1:
                        continue
                    if ac not in ac2isoformlist:
                        ac2isoformlist[ac] = []
                    ac2isoformlist[ac].append(isoform)
                prd = "<http://glygen-vm-prd.biochemistry.gwu.edu/ontology/canonical>"
                if line.find(prd) != -1:
                    isoform = line.strip().split(" ")[0].split("/")[-1][:-1]
                    flag = line.strip().split(" ")[-2].strip().split("^")[-0].replace("\"", "")
                    if flag == "true":
                        is_canon[isoform] = True
       
        newcanon2isoformlist = {}
        for ac in ac2isoformlist:
            for isoform in ac2isoformlist[ac]:
                if isoform in is_canon:
                    canon = isoform
                    newcanon2isoformlist[canon] = ac2isoformlist[ac]
        #for canon in newcanon2isoformlist:
        #    print canon, sorted(set(newcanon2isoformlist[canon]))
        #sys.exit()



	client = MongoClient('mongodb://localhost:27017')
	db = client["glygen"]
        query_obj = {"species.taxid": {'$eq': int(tax_id)}}

        for doc in db["c_protein"].find(query_obj):
            canon = doc["uniprot_canonical_ac"]
            isoform_list = []
            for obj in doc["isoforms"]:
                isoform_list.append(obj["isoform_ac"])
            if canon in newcanon2isoformlist:
                if sorted(set(newcanon2isoformlist[canon])) != sorted(set(isoform_list)):
                    print canon



if __name__ == '__main__':
	main()


