import os,sys
import json
import csv
import time

import commands
import glob

from optparse import OptionParser

import libgly




def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog ")
    parser.add_option("-v","--ver",action="store",dest="ver",help="data release version")

    (options,args) = parser.parse_args()
    for file in ([options.ver]):
        if not (file):
            parser.print_help()
            sys.exit(0)

   
    ver = options.ver
    #rel_dir = "/data/shared/glygen/releases/data/v-%s/" % (ver)
    rel_dir = "."

    file_list_new = glob.glob("unreviewed/*")
    file_list_old = glob.glob(rel_dir + "/reviewed/*")


    seen_old, seen_new = {}, {}
    for in_file in file_list_new:
        file_name = in_file.split("/")[-1]
        seen_new[file_name] = True
    for in_file in file_list_old:
        file_name = in_file.split("/")[-1]
        if file_name.find("GLY_") == -1 and file_name.find(".stat.csv") == -1:
            #file_name = file_name.replace("citations_glycosylation_", "citations_glycosylation_sites_")
            #file_name = file_name.replace("citations_phosphorylation_", "citations_phosphorylation_sites_")
            #file_name = file_name.replace("citations_glycation_", "citations_glycation_sites_")
            file_name = file_name.replace("_biomarkers_cancer", "_biomarkers")
            file_name = file_name.replace("_o_glcnac_", "_oglcnac_")
            seen_old[file_name] = True


    for file_name in list(set(seen_new.keys() + seen_old.keys())):
        print file_name, file_name in seen_new, file_name in seen_old

 
           
 
if __name__ == '__main__':
            
    main()



