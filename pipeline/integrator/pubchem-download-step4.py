import os
import csv
import sys
import json
import commands
import glob
import gzip
from optparse import OptionParser

sys.path.append('../../glytools/')
import libgly



def main():


    cmd = "cat generated/pubchem/compound/header.csv "
    cmd += " generated/pubchem/compound/cid2inchi.batch.*.csv "
    cmd += " > generated/pubchem/compound/cid2inchi.csv"
    x = commands.getoutput(cmd)
    #print cmd

                


if __name__ == '__main__':
        main()
