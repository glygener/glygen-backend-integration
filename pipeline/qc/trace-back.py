#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
import random



__version__="1.0"
__status__ = "Dev"


#######################################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
	root_obj =  config_obj[config_obj["server"]]["rootinfo"]

	db_obj = {}
	taxid2name = {9606:"Homo sapiens", 10090:"Mus musculus"}
	
        dst_dir = "/data/projects/glygen/generated/datasets/reviewed/"
        src_dir = "/data/projects/glygen/downloads/uckb/"

        paired_files = [
                { 
                    "src": src_dir + "/human_query_results_2018_08_07_19_13_12.txt"
                    ,"dst": dst_dir + "/human_proteoform_glycosylation_sites_unicarbkb_glytoucan.csv"
                }
        ]

        src_grid = []
        with open(paired_files[0]["src"], 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                row_count += 1
                if row_count == 1:
                    field_list = row
                else:
                    row_obj = {}
                    for j in xrange(0,len(field_list)):
                        field_name  = field_list[j]
                        row_obj[field_name] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                    src_grid.append(row_obj)
                
        for row_obj in src_grid:
            print row_obj





if __name__ == '__main__':
        main()



