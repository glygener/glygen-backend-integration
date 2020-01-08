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

    in_file = "/data/projects/glygen/generated/misc/dataset-masterlist.json"
    ds_obj_list = json.loads(open(in_file, "r").read())
    ntested, npassed, nfailed = 0, 0, 0
    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if ds_format != "csv":
            continue
        file_list = []
        for species in obj["categories"]["species"]:
            f = "unreviewed/%s_%s_%s.%s" % (species, mol, ds_name, ds_format)
            file_list.append(f)
        if file_list == []:
            f = "unreviewed/%s_%s.%s" % (mol, ds_name, ds_format)
            file_list.append(f)

        for in_file in file_list:
            if os.path.isfile(in_file) == True:
                data_frame = {}
                libgly.load_sheet(data_frame, in_file, ",")
                f_list = data_frame["fields"]
                n_fields = len(f_list)
                flag = True
                for row in data_frame["data"]:
                    n_cols = len(row)
                    if n_fields != n_cols:
                        flag = False
                        print "Bad row"
                        print in_file
                        print f_list
                        print row
                        break
                ntested += 1
                npassed += 1 if flag == True else 0
                nfailed += 1 if flag == False else 0
                print "%s tested, %s passed, %s failed" % (ntested,npassed, nfailed)
















if __name__ == '__main__':
        main()


