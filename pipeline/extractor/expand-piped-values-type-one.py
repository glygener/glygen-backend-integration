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





def expand_df(df, expand_j):
    
    new_df = []
    for row in df:
        for expand_val in row[expand_j].split("|"):
            new_row = []
            for j in xrange(0, len(row)):
                val = expand_val if j == expand_j else row[j]
                new_row.append(val)
            new_df.append(new_row)

    return new_df


#######################################
def main():


    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " )
    parser.add_option("-i","--infile",action="store",dest="infile",help="Input file")
    parser.add_option("-f","--fieldlist",action="store",dest="fieldlist",help="field_a,field_b")

    (options,args) = parser.parse_args()
    for file in ([options.infile, options.fieldlist]):
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
    field_list = options.fieldlist.split(",")


    sheet_obj = {}
    libgly.load_sheet(sheet_obj, in_file, ",")

    df = sheet_obj["data"]

    print "\"%s\"" % ("\",\"".join(sheet_obj["fields"]))
    for field in field_list:
        field_ind = sheet_obj["fields"].index(field)
        df = expand_df(df, field_ind)
    
    for row in df:
        print "\"%s\"" % ("\",\"".join(row))






if __name__ == '__main__':
        main()



