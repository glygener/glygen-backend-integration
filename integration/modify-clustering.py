import os,sys
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


######################
def load_dataframe(data_frame, sheet_name, in_file, separator):

    data_frame[sheet_name] = {}
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='"')
        row_count = 0
        for row in csv_grid:
            row_count += 1
            if row_count == 1:
                field_list = row
            else:
                row_obj = {}
                for j in xrange(1,len(field_list)):
                    field_name  = field_list[j]
                    row_obj[field_name] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                main_id = row[0].strip()
                if main_id not in data_frame[sheet_name]:
                    data_frame[sheet_name][main_id] = []
                data_frame[sheet_name][main_id].append(row_obj)
    return




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

    gene_file = "downloads/protein/UP000005640_9606_acc2groupid_2017_11.tsv"
    if tax_id == "10090":
        gene_file = "downloads/protein/UP000000589_10090_acc2groupid_2017_11.tsv"

    ac2genename = {}
    seen = {"genename":{}} 
    with open(gene_file, "r") as FR:
        for line in FR:
            ac = line.strip().split("\t")[0].split("-")[0]
            gene_name = line.strip().split("\t")[-1].lower()
            seen["genename"][gene_name] = True
            ac2genename[ac] = line.strip().split("\t")[-1]



    data_frame = {}

    in_file = "reviewed/%s_protein_idmapping.csv" % (species)
    load_dataframe(data_frame, "idmapping", in_file, ",")
    in_file = "reviewed/%s_protein_information.csv" % (species)
    load_dataframe(data_frame, "information", in_file, ",")
    

    #check sanity of gene names
    for canon in data_frame["idmapping"]:
        for obj in data_frame["idmapping"][canon]:
            for gene_name in obj["gene_name"]:
                #if 0 and gene_name.lower() not in seen["genename"]:
                if gene_name.lower() not in seen["genename"]:
                    ac = canon.split("-")[0]
                    #print gene_name, canon, ac2genename[ac]
                    obj["gene_name"] = [ac2genename[ac]]



    canon2genename, canon2status = {}, {}
    cls_dict = {"reviewed":{}, "unreviewed":{}}
    na_count = 1
    gene2canon = {}
    for canon in data_frame["idmapping"]:
        cls_dict["reviewed"][canon], cls_dict["unreviewed"][canon] = [], []
        for obj in data_frame["idmapping"][canon]:
            cls_dict["reviewed"][canon] += obj["reviewed_isoforms"] 
            cls_dict["unreviewed"][canon] += obj["unreviewed_isoforms"]
            for gene_name in obj["gene_name"]:
                if gene_name.strip() == "-":
                    gene_name = "UNKNOWN_%s" % (na_count)
                    na_count += 1
                if gene_name not in gene2canon:
                    gene2canon[gene_name] = []
                gene2canon[gene_name].append(canon)
                canon2genename[canon] = gene_name
            canon2status[canon] = obj["status"][0]
                                    
     
    canon2remove = {}
    canons4merge = {}
    for gene_name in gene2canon:
        if len(gene2canon[gene_name]) > 1:
            longest_canon = ""
            max_len = 0
            for canon in gene2canon[gene_name]:
                #print gene_name, canon, canon2status[canon]
                canon2remove[canon] = True
                canon_len = data_frame["information"][canon][0]["protein_length"][0]
                if canon_len > max_len:
                    longest_canon = canon
                    max_len = int(canon_len)
            canons4merge[longest_canon] = gene2canon[gene_name]
            canons4merge[longest_canon].remove(longest_canon)
            #print longest_canon, canons4merge[longest_canon]


    for canon in data_frame["idmapping"]:
        for obj in data_frame["idmapping"][canon]:
            row = [
                canon
                ,"|".join(obj["status"])
                ,"|".join(obj["gene_name"])
                ,"|".join(obj["reviewed_isoforms"])
                ,"|".join(obj["unreviewed_isoforms"])

            ]
            if canon not in canon2remove:
                print "%s,%s,%s,%s,%s"  % (row[0],row[1],canon2genename[canon],row[3], row[4])
   

    for longest_canon in canons4merge:
        for canon in canons4merge[longest_canon]:
            #print "AAA",longest_canon,canon,canon2status[longest_canon],canon2status[canon]
            cls_dict["reviewed"][longest_canon] += cls_dict["reviewed"][canon]
            cls_dict["unreviewed"][longest_canon] += cls_dict["unreviewed"][canon]
        row = [
                longest_canon
                ,canon2status[longest_canon]
                ,canon2genename[longest_canon]
                ,"|".join(sorted(set(cls_dict["reviewed"][longest_canon])))
                ,"|".join(sorted(set(cls_dict["unreviewed"][longest_canon])))
       ]
        print "%s,%s,%s,%s,%s"  % (row[0],row[1],row[2],row[3], row[4])





if __name__ == '__main__':
	main()
