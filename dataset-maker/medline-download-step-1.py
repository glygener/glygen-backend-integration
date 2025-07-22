import os
import csv
import sys
import json
import glob
import gzip
import xml.etree.ElementTree as ET
from Bio import SeqIO
import time
import requests


import libgly3





def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
 
    global path_obj
    global downloaded_list

    path_obj = config_obj["pathinfo"]

    DEBUG = False
    #DEBUG = True


    file_list = glob.glob(path_obj["downloads"] + "ncbi/medline_txt/pmid.*.txt")
    file_list += glob.glob(path_obj["downloads"] + "ncbi/medline_json/pmid.*.json")

    is_downloaded = {}
    for in_file in file_list:
        pmid = in_file.split(".")[-2]
        is_downloaded[pmid] = True


    FW = open("logs/medline-download-list.txt", "w")
    seen = {}
    file_list = glob.glob("unreviewed/*.csv")
    file_list += glob.glob("unreviewed/*.tsv") 

    for in_file in file_list:
        data_frame = {}
        libgly3.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        if "xref_key" not in f_list:
            continue
        for row in data_frame["data"]:
            xref_key, xref_id = row[f_list.index("xref_key")], row[f_list.index("xref_id")]
            if xref_key.find("pubmed") == -1:
                continue
            if xref_id == "" or xref_id.isdigit() == False:
                continue
            if xref_id not in seen and xref_id not in is_downloaded:        
                FW.write("%s\n" % (xref_id))
                seen[xref_id] = True
    FW.close()  

            
    return


                


if __name__ == '__main__':
        main()
