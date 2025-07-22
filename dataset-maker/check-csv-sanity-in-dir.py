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


    file_list = glob.glob("unreviewed/*.csv")
    #file_list = glob.glob("tmp/*.csv")
    for in_file in file_list:
        if os.path.isfile(in_file) == True:
            with open(in_file, "r") as FR:
                lcount = 0
                n_fields = 0
                f_list = []
                pass_flag = True
                failed_rows = 0
                for line in FR:
                    lcount += 1
                    row = line.strip().split("\",\"")
                    if lcount == 1:
                        n_fields = len(row)
                        f_list = row
                    else:
                        n_values = len(row)
                        if n_fields != n_values:
                            failed_rows += 1
                            print "Bad row in %s, row number=%s" % (in_file,lcount)
                            print "ncols=%s, nvalues=%s" % (n_fields, n_values)
                            exit()

                if failed_rows == 0:
                    print "passed ", in_file













if __name__ == '__main__':
        main()


