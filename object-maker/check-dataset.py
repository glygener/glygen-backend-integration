#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import subprocess
from optparse import OptionParser

import csvutil



def sort_release_list(tmp_list, reversed_flag):

    factor_list = [100000000, 1000, 1]
    rel_dict = {}
    for rel in tmp_list:
        parts = rel.split(".")
        ordr = 0
        for i in range(0,len(parts)):
            ordr += factor_list[i]*int(parts[i])
        rel_dict[ordr] = rel
    
    release_list = []

    for ordr in sorted(rel_dict, reverse=reversed_flag):
        release_list.append(rel_dict[ordr])

    return release_list


def get_record_count_detailed(in_file):
    
    file_ext = in_file.split(".")[-1].lower()
    field_count, row_count, id_count = 1, 1, 1
    
    res_obj = {}
    if file_ext in ["csv", "tsv"]:
        sep = "\t" if file_ext == "tsv" else ","
        res_obj = csvutil.get_sheet_stats_detailed(in_file, sep)
    elif file_ext in ["fasta"]:
        res_obj["rowcount"] = len(list(SeqIO.parse(in_file, "fasta")))

    return res_obj


def get_record_count(in_file):

    file_ext = in_file.split(".")[-1].lower()
    field_count, row_count, id_count = 1, 1, 1
    if file_ext in ["csv", "tsv"]:
        sep = "\t" if file_ext == "tsv" else ","
        field_count, row_count, id_count = csvutil.get_sheet_stats(in_file, sep)
    elif file_ext in ["fasta"]:
        id_count = len(list(SeqIO.parse(in_file, "fasta")))

    return field_count, row_count, id_count


def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version=" ")
    parser.add_option("-i","--infile",action="store",dest="infile",help="Unreviewed dataset file")
    parser.add_option("-r","--refrel",action="store",dest="refrel",help="reference release")
    (options,args) = parser.parse_args()
    for file in ([options.infile, options.refrel]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    global data_release_dir
    global unreviewed_dir
    global current_rel

    current_rel = "x.x.x"
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    data_release_dir = "/data/shared/glygen/releases/data/"
    unreviewed_dir = wrk_dir + "/unreviewed/"

    in_file = options.infile
    ref_rel = options.refrel
    
    file_name = in_file.split("/")[-1]
    old_in_file = data_release_dir + "v-" + ref_rel + "/reviewed/" + file_name

    #field_count, row_count, id_count = get_record_count(old_in_file)
    #print (field_count, row_count, id_count)
    #field_count, row_count, id_count = get_record_count(in_file)
    #print (field_count, row_count, id_count)

    res_obj_one = get_record_count_detailed(old_in_file)
    res_obj_two = get_record_count_detailed(in_file)

    f_list_one = list(res_obj_one["uniquecount"].keys())
    f_list_two = list(res_obj_two["uniquecount"].keys())

    row = ["rowcount", res_obj_one["rowcount"], res_obj_two["rowcount"]]
    print (row)
    for f in sorted(set(f_list_one + f_list_two)):
        v_one = res_obj_one["uniquecount"][f] if f in res_obj_one["uniquecount"] else ""
        v_two = res_obj_two["uniquecount"][f] if f in res_obj_two["uniquecount"] else ""
        row = [f, v_one, v_two]
        print (row)


if __name__ == '__main__':
    main()

