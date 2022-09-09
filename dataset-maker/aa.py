import os,sys
import json
import csv

from optparse import OptionParser

import commands
import glob
import csvutil


__version__="1.0"
__status__ = "Dev"



###############################
def main():


    in_file = ""
    data_frame = {}
    csvutil.load_sheet(data_frame, in_file, [], ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:


















if __name__ == '__main__':
        main()


