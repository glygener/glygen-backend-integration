import os
import sys
import gzip
import json
from optparse import OptionParser
import glob
import pubmed_parser as pp






###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " )
    parser.add_option("-s","--start",action="store",dest="start",help="")
    parser.add_option("-e","--end",action="store",dest="end",help="")

    (options,args) = parser.parse_args()
    for file in ([options.start, options.end]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    start = int(options.start)
    end = int(options.end)

    xml_folder = "downloads/ncbi/medline_xml/"
    out_folder = "downloads/ncbi/medline_list/"
    gz_xml_file_list = json.load(open("logs/medline_xml_file_list.json"))
    gz_xml_file_list = sorted(gz_xml_file_list)

    for gz_xml_file in gz_xml_file_list[start-1:end]:
        file_name = gz_xml_file.split("/")[-1].replace(".xml.gz", "")
        out_file = out_folder + file_name + ".txt"
        if os.path.isfile(out_file) == True:
            continue
        dicts_out = pp.parse_medline_xml(gz_xml_file,
            year_info_only=False,
            nlm_category=False,
            author_list=True,
            reference_list=False
        )
        tmp_list = []
        for obj in dicts_out:
            tmp_list.append(obj["pmid"])
        with open(out_file, "w") as FW:
            FW.write("%s\n" % ("\n".join(tmp_list))) 


if __name__ == '__main__':
    main()

