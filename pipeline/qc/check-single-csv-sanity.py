import os,sys
import json
import csv

from optparse import OptionParser

import commands
import glob

sys.path.append('../../glytools/')
import libgly


__version__="1.0"
__status__ = "Dev"



###############################
def main():

    in_file = "downloads/unicarbkb/human29112019.csv"
    if os.path.isfile(in_file) == True:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        n_fields = len(f_list)
        flag = True
        row_count = 0
        for row in data_frame["data"]:
            row_count += 1
            n_cols = len(row)
            if n_fields != n_cols:
                flag = False
                print "Bad row, row number=%s" % (row_count)
                print in_file
                print f_list
                print row
                break
                














if __name__ == '__main__':
        main()


