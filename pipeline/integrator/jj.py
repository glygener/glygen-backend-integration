import os
import csv
import sys
import json
import commands
import glob
import gzip
import xml.etree.ElementTree as ET
from Bio import SeqIO
import time
import requests


sys.path.append('../../glytools/')
import libgly

def extract_pmid_list(file_list, field_name, delim, pmid_list, downloaded_list):

    tmp_list = []
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, delim)
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            if field_name == "xref_id":
                if row[f_list.index("xref_key")].find("pubmed") == -1:
                    continue
            val = row[f_list.index(field_name)].replace(",", " ").replace("  ", " ")
            val = val.replace(";", " ").replace("  ", " ")
            id_list = val.split(" ")
            for pmid in id_list:
                pmid = pmid.strip()
                if pmid == "" or pmid.isdigit() == False:
                    continue
                tmp_list.append(pmid)

        set_one = set(tmp_list) - set(pmid_list)
        set_two = set_one - set(downloaded_list)
        pmid_list += list(set_two)
        n1 = len(tmp_list)
        n2 = len(set_two)



    return





def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
 
    global path_obj
    global downloaded_list
    global pmid_list

    path_obj = config_obj["pathinfo"]

    script_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    downloaded_list = []
    pmid_list = ["PMC6243375"]
    for pmid in pmid_list:
        time.sleep(1)
        out_file = "pmid.%s.txt" % (pmid)
        url = "%s?db=pmc&rettype=medline&id=%s"  % (script_url, pmid)
        res = requests.get(url)
        if res.content.find("PMID") != -1:
            with open(out_file, 'w') as FW:
                FW.write("%s\n" % (res.content))


    


            
    return


                


if __name__ == '__main__':
        main()
