import os
import sys
import gzip
import json
from optparse import OptionParser
import glob
import pubmed_parser as pp
import subprocess



def load_pmid_dict(in_file):

    seen = {}
    pmid_list =  open(in_file, "r").read().split("\n")
    for pmid in pmid_list:
        if pmid.strip() == "":
            continue
        seen[pmid] = True

    return seen


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

    json_folder = "downloads/ncbi/medline_json/"
    list_file = "logs/medline-download-list.txt"
    pmid_dict = load_pmid_dict(list_file)


    count = 0
    logfile = "logs/xml2json.%s-%s.log" % (start, end)
    with open(logfile, "w") as FL:
        FL.write("Started processing\n")

    xml_file_list = json.loads(open("logs/xml_file_list.json", "r").read())
    xml_file_list = sorted(xml_file_list)
    for xml_file in xml_file_list[start-1:end]:
        dicts_out = pp.parse_medline_xml(xml_file,
            year_info_only=False,
            nlm_category=False,
            author_list=True,
            reference_list=False
        )
        doc = {"TI":"", "JT":"", "BTI":"", "DP":"", "AU":""}
        found = 0
        for obj in dicts_out:
            doc = {"PMID":obj["pmid"], "TI":"", "JT":"", "BTI":"", "DP":"", "AU":""}
            doc["TI"] = obj["title"]
            doc["DP"] = obj["pubdate"]
            doc["JT"] = obj["journal"]
            tmp_list = []
            for o in obj["authors"]:
                f = o["forename"] if "forename" in o else ""
                l = o["lastname"] if "lastname" in o else ""
                tmp_list.append(l + " " + f)
            doc["AU"] = ", ".join(tmp_list)
            pmid = doc["PMID"]
            if pmid.strip() != "" and pmid in pmid_dict:
                out_file = json_folder + "pmid.%s.json" % (pmid)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(doc, indent=4)))
                found += 1 
        with open(logfile, "a") as FL: 
            FL.write("finished %s, found %s PMIDs\n" % (xml_file, found))


if __name__ == '__main__':
    main()

