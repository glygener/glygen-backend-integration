import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient

sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"



###############################
def main():

    pattern = "unreviewed/*.csv"
    for in_file in glob.glob(pattern):
        seen = {}
        out_buffer = ""
        flag = False
        with open(in_file, "r") as FR:
            for line in FR:
                if line.strip() not in seen:
                    out_buffer += line
                else:
                    flag = True
                seen[line.strip()] = True
        if flag == True:
            out_file = in_file
            with open(out_file, "w") as FW:
                FW.write("%s" % (out_buffer))
            print "Cleaned %s " % (in_file)





if __name__ == '__main__':
	main()

