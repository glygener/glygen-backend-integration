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

    in_file = "generated/misc/field_names.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    seen = {}
    for row in data_frame["data"]:
        seen[row[0]] = True
    

    config_obj = json.loads(open("/data/projects/glygen/generated/misc/dataset-masterlist.json", "r").read())
   
    for species in ["human", "mouse", "rat"]:
        for cat in ["protein", "glycan", "proteoform"]:
            ds_list = config_obj[cat]["common"]
            ds_list += config_obj[cat][species] if species in config_obj[cat] else []

            for ds in ds_list:
                if ds in ["allsequences", "canonicalsequences"]:
                    continue
                ext = "csv"
                file_name = "%s_%s_%s.%s" % (species, cat, ds, ext)
                in_file = "unreviewed/%s" % (file_name)
                if os.path.isfile(in_file) == True:
                    data_frame = {}
                    libgly.load_sheet(data_frame, in_file, ",")
                    f_list = data_frame["fields"]
                    for f in f_list:
                        if f not in seen:
                            print "%s,%s" % (f, file_name)

















if __name__ == '__main__':
        main()


