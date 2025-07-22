import os
import sys
import gzip
import json
import glob
import pubmed_parser as pp
import re

import spacy
from optparse import OptionParser




def load_pmid_dict(in_file):

    seen = {}
    pmid_list =  open(in_file, "r").read().split("\n")
    for pmid in pmid_list:
        if pmid.strip() == "":
            continue
        seen[pmid] = True

    return seen




def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-s","--start",action="store",dest="start",help="")
    parser.add_option("-e","--end",action="store",dest="end",help="")

    (options,args) = parser.parse_args()
    for key in ([options.start, options.end]):
        if not (key):
            parser.print_help()
            sys.exit(0)


    glygen_dict = load_pmid_dict("logs/medline-glygen-list.txt")


    start_idx = int(options.start)
    end_idx = int(options.end)
    
    xml_dir = "downloads/ncbi/medline_xml/"
    out_dir = "downloads/ncbi/medline_abstracts/"
    log_file = "logs/extract_abstracts.%s.%s.log" % (start_idx, end_idx)

    with open(log_file, "w") as FL:
        FL.write("Started logging\n")

    file_list = sorted(glob.glob(xml_dir + "/*.xml.gz"))
    end_idx = len(file_list) if end_idx > len(file_list) else end_idx

    for gz_xml_file in file_list[start_idx-1:end_idx]:
        if not gz_xml_file.endswith('.gz'):
            continue
        file_name = gz_xml_file.split("/")[-1].replace(".xml.gz", "")
        dicts_out = pp.parse_medline_xml(gz_xml_file)
        n_found = 0
        idx = 0
        for obj in dicts_out:
            if "pmid" not in obj:
                continue
            if obj["pmid"] == "":
                continue
            if "abstract" not in obj:
                continue
            if obj["abstract"] == "":
                continue
            doc_id = obj["pmid"]
            if doc_id not in glygen_dict:
                continue            

            if idx > 0 and idx%1000 == 0:
                with open(log_file, "a") as FL:
                    FL.write("checked %s PMIDs, %s found from %s\n" % (idx, n_found, gz_xml_file))
            idx += 1
            title_text, abstract_text = obj["title"], obj["abstract"]
            mesh = obj["mesh_terms"]
            publication_types = obj["publication_types"]
            text = "{} {}".format(title_text, abstract_text).strip()
            doc = {"docId":doc_id, "text":text} 
            if len(title_text) > 0:
                doc["title"] = {"charStart": 0, "charEnd": len(title_text) - 1}
            if len(publication_types) > 0:
                doc["type"] = publication_types
            if len(mesh) > 0:
                doc["mesh"] = mesh
            out_file = out_dir + "pmid.%s.txt" % (doc_id)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (doc["text"]))
        with open(log_file, "a") as FL:
            FL.write("final %s found from %s\n" % (n_found, gz_xml_file))
        
    return




if __name__ == '__main__':
    main()

