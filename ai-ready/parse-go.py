import os
import sys
import json
import glob
from optparse import OptionParser


import libgly



def load_goid2lineage(in_file):
    
    goid2lineage = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        go_id, go_ns, lineage = row[0], row[1], row[2].split(";")
        
        if go_id not in goid2lineage:
            goid2lineage[go_id] = {}
        goid2lineage[go_id][go_ns] = lineage

    return goid2lineage

def load_canon2goid(in_file, lineage_file, depth):

    goid2lineage = load_goid2lineage(lineage_file)

    canon2goid = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon, go_id = row[0], row[2].replace("_", ":")
        if go_id not in goid2lineage:
            continue
        for go_ns in goid2lineage[go_id]:
            lineage = goid2lineage[go_id][go_ns]
            if len(lineage) > depth:
                last_p = lineage[-depth]
                canon2goid[canon] = last_p + "|" + go_ns

    return canon2goid


def main():

    global go_dict

    lineage_file = "outdir/go_lineage.csv"
    ann_file = "unreviewed/human_protein_go_annotation.csv"
    depth = 3
    canon2goid = load_canon2goid(ann_file, lineage_file, depth)
    
    for canon in canon2goid:
        print canon, canon2goid[canon]

    
    return


                


if __name__ == '__main__':
        main()
