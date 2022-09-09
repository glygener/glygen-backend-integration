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

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse/rat") 
    parser.add_option("-r","--release",action="store",dest="release",help="1.0.13")


    (options,args) = parser.parse_args()
    for file in ([options.species, options.release]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    species = options.species
    release = options.release

    config_obj = json.loads(open("/data/projects/glygen/generated/misc/dataset-masterlist.json", "r").read())
    
    for cat in ["protein", "glycan", "proteoform"]:
        ds_list = config_obj[cat]["common"]
        ds_list += config_obj[cat][species] if species in config_obj[cat] else []


        for ds in ds_list:
            if ds in ["allsequences", "canonicalsequences"]:
                continue
            ext = "csv"
            file_name = "%s_%s_%s.%s" % (species, cat, ds, ext)
            in_file_one = "unreviewed/%s" % (file_name)
            in_file_two = "/data/projects/glygen/releases/data/v-%s/reviewed/%s" % (release,file_name)
            if os.path.isfile(in_file_one) == False or os.path.isfile(in_file_two) == False:
                print "%s was not found in both versions" % (file_name) 
                continue
            data_frame_one, data_frame_two = {},  {}
            libgly.load_sheet(data_frame_one, in_file_one, ",")
            libgly.load_sheet(data_frame_two, in_file_two, ",")
            f_list_one = data_frame_one["fields"]
            f_list_two = data_frame_two["fields"]
            if f_list_one != f_list_two:
                set_one = set(f_list_one)
                set_two = set(f_list_two)
                diff_one = "|".join(list(set_one - set_two))
                diff_two = "|".join(list(set_two - set_one))
                print "%s,%s,%s" %(file_name, diff_one, diff_two)

















if __name__ == '__main__':
        main()


