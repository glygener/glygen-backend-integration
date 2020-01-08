import os,sys
import json
import csv

from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


import commands
import glob
import sparqlutil
#import requests
import gzip


sys.path.append('../../glytools/')
import libgly


__version__="1.0"
__status__ = "Dev"




def make_taxid2name_ds():

    seen = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/taxa.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        tax_id = row[f_list.index("TaxID")]
        seen[tax_id] = True

    newrow = ["tax_id", "tax_name"]
    print "\"%s\"" % ("\",\"".join(newrow))

    in_file = path_obj["downloads"] +  "/ncbi/taxonomy/names.dmp"
    with open(in_file, "r") as FR:
        for line in FR:
            parts = line.strip().split("|")
            if parts[3].strip() == "scientific name":
                tax_id = parts[0].strip()
                tax_name = parts[1].strip()
                if tax_id in seen:
                    newrow = [tax_id, tax_name]
                    print "\"%s\"" % ("\",\"".join(newrow))
    return



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-d","--dataset",action="store",dest="dataset",help="[masterlist, transcriptlocus]")

    (options,args) = parser.parse_args()
    for file in ([options.dataset]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    global config_obj
    global sparql
    global graph_uri
    global prefixes
    global data_grid
    global species_obj


    dataset = options.dataset

    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]

    global path_obj
    path_obj = config_obj["pathinfo"]


    data_grid = {}
    if dataset == "taxid2name":
        make_taxid2name_ds()






if __name__ == '__main__':
        main()


