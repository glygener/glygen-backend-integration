import os
import sys
import json
import glob
from optparse import OptionParser


import libgly


def main():

    in_file = "tmp/toy.csv"
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    row = ["sample","cls"]
    print "\"%s\"" % ("\",\"".join(row))
    for f in f_list[1:]:
        row = [f,f]
        print "\"%s\"" % ("\",\"".join(row))



            
    return


                


if __name__ == '__main__':
        main()
