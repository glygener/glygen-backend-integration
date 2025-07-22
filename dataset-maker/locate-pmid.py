import os
import sys
import gzip
import json
from optparse import OptionParser
import glob
import pubmed_parser as pp
import subprocess



def load_pmid_list(in_file):

    seen = {}
    pmid_list =  open(in_file, "r").read().split("\n")
    for pmid in pmid_list:
        if pmid.strip() == "":
            continue
        seen[pmid] = True

    return list(seen.keys())


###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")

    parser.add_option("-i","--infile",action="store",dest="infile",help="")
    parser.add_option("-s","--start",action="store",dest="start",help="")
    parser.add_option("-e","--end",action="store",dest="end",help="")

    (options,args) = parser.parse_args()
    for key in ([options.infile, options.start, options.end]):
        if not (key):
            parser.print_help()
            sys.exit(0)
            
    global data_dir

    list_file = options.infile
    start = int(options.start)
    end = int(options.end)

    data_dir = "/data/shared/nlp/"
    pmid_list = load_pmid_list(list_file)

    xml_folder = data_dir + "medline/xml/"
    res_folder = data_dir + "medline/json/"
    pmid_list = open(list_file, "r").read().split("\n")

    count = 0
    logfile = "logs/locate.%s.%s.log" % (start, end)
    with open(logfile, "w") as FL:
        FL.write("Started processing\n")

    file_list = sorted(glob.glob(xml_folder + "/*.xml.gz"))
    found = 0
    for idx in range(start-1, end):
        gz_xml_file = file_list[idx]
        if not gz_xml_file.endswith('.gz'):
            continue
        for pmid in pmid_list:
            if pmid.strip() == "":
                continue
            cmd = "gunzip -c " + gz_xml_file + " | grep " + pmid
            x = subprocess.getoutput(cmd)
            if x.strip() != "":
                with open(logfile, "a") as FL: 
                    FL.write(" ... %s found in %s\n" % (pmid, gz_xml_file))
        with open(logfile, "a") as FL:
            FL.write(" ... finished parsing %s\n" % (gz_xml_file))


if __name__ == '__main__':
    main()

