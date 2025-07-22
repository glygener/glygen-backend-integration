import os
import sys
import gzip
import json
from optparse import OptionParser
import glob
import subprocess
import csv
import io
import libgly



###############################
def main():




    seen = {}
    for species in ["mouse", "rat", "fruitfly", "yeast"]:
        in_file = "downloads/alliance_genome/current/%s_disease_genome_alliance.tsv" % (species)
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, "\t")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            if len(row) != len(f_list):
                continue
            do_id = row[f_list.index("DOID")].split(":")[1]
            seen[do_id] = "A"

    in_file = "unreviewed/protein_disease_idmap.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")].strip()
        if do_id not in seen:
            seen[do_id] = "B"
        else:
            seen[do_id] += "B"


    for do_id in seen:
        print seen[do_id], do_id




if __name__ == '__main__':
    main()

