#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict


sys.path.append('../../glytools/')
import libgly







#######################################
def main():


    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " )
    parser.add_option("-i","--infile",action="store",dest="infile",help="Input file")

    (options,args) = parser.parse_args()
    for file in ([options.infile]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    global config_obj
    global sparql
    global graph_uri
    global prefixes
    global data_grid

    config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())


    in_file = options.infile


    sheet_obj = {}
    libgly.load_sheet(sheet_obj, in_file, ",")


    print "\"%s\"" % ("\",\"".join(sheet_obj["fields"]))
        
    for row in sheet_obj["data"]:
        df = []
        max_n = 0
        for val in row:
            val_list = val.split("|")
            df.append(val_list)
            if len(val_list) > max_n:
                max_n = len(val_list)
        new_df = []
        for i in xrange(0, max_n):
            new_df.append([])
            for j in xrange(0, len(row)):
                val = df[j][i] if i < len(df[j]) else new_df[i-1][j]
                new_df[i].append(val)
        for row in new_df:
            print "\"%s\"" % ("\",\"".join(row))







if __name__ == '__main__':
        main()



