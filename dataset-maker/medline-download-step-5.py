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

    xml_folder = "downloads/ncbi/medline_xml/"
    list_folder = "downloads/ncbi/medline_list/"
    json_folder = "downloads/ncbi/medline_json/"
    list_file = "logs/medline-download-list.txt"
    pmid_dict = load_pmid_dict(list_file)



    xml_file_list = json.loads(open("logs/xml_file_list.json", "r").read())
    xml_file_list = sorted(xml_file_list)
    
    n = len(xml_file_list)
    batch_size = int(len(xml_file_list)/10) if len(xml_file_list) > 10 else len(xml_file_list)

   
    range_list = [] 
    for i in range(0, 12):
        s = i*batch_size + 1
        e = s + batch_size
        if e >= n:
            e = n
            range_list.append({"s":s, "e":e})
            break
        range_list.append({"s":s, "e":e}) 

    for o in range_list:
        cmd = "nohup python3 extract-medline.py -s %s -e %s  & " % (o["s"], o["e"])
        os.system(cmd)



if __name__ == '__main__':
    main()

