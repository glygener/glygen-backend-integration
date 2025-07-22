#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import subprocess
from optparse import OptionParser


def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version=" ")
    parser.add_option("-r","--refrel",action="store",dest="refrel",help="reference release")
    (options,args) = parser.parse_args()
    for file in ([options.refrel]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    ref_rel = options.refrel
    d_one = "/data/shared/glygen/releases/data/v-%s/jsondb/" % (ref_rel)
    d_two = "jsondb/"

    seen = {}
    d_list = glob.glob(d_one + "*")
    for d in glob.glob(d_one + "*")  + glob.glob(d_two + "*"):
        dd = d.split("/")[-1]
        if dd[-2:] == "db":
            seen[dd] = True
    
    for dd in seen:
        n1, n2 = 0, 0
        dd_one = d_one + dd 
        dd_two = d_two + dd
        if os.path.isdir(dd_one):
            n1 = len(glob.glob(dd_one + "/*.json"))
        if os.path.isdir(dd_two):
            n2 = len(glob.glob(dd_two + "/*.json"))
        print (dd, n1, n2)




            

if __name__ == '__main__':
    main()

