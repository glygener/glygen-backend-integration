#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import subprocess
from optparse import OptionParser

import csvutil




def main():

    for record_type in ["protein", "glycan"]:
        expected_dslist =  csvutil.get_expected_dslist(record_type, [], {})
        mising_files = []
        for file_name in expected_dslist:
            file_ext = "fasta" if file_name.find("_allsequences") != -1 else "csv"
            in_file = "reviewed/%s.%s" % (file_name,file_ext)
            if os.path.isfile(in_file) == False:
                mising_files.append(in_file)
                print (file_name, "missing dataset")

    if mising_files != []:
        print ("The following files are missing:")
        print (json.dumps(mising_files, indent=4))
        exit()
    else:
        print ("All required dataset files are available")


if __name__ == '__main__':
    main()

