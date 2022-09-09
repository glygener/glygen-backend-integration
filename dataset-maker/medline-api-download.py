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
            val = ""
            if field_name in f_list:
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



def extract_pmid_list_from_xref(pmid_list, downloaded_list):

    tmp_list  = []
    for gpff_file in glob.glob("downloads/ncbi/refseq/current/refseq_protein_all_7227.gpff"):
        for record in SeqIO.parse(gpff_file, "genbank"):
            if "references" not in record.annotations:
                continue
            for ref in record.annotations["references"]:
                pmid = ref.pubmed_id
                if pmid != "":
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


    downloaded_list = []
    pmid_list = []

    file_list = glob.glob(path_obj["downloads"] + "ncbi/medline/pmid.*.txt")
    for in_file in file_list:
        pmid = in_file.split(".")[-2]
        downloaded_list.append(pmid)
    



    #Glytoucan publications
    file_list = [path_obj["downloads"] +  "glytoucan/current/export/pubs.tsv"]
    extract_pmid_list(file_list, "PubMedID", "\t", pmid_list, downloaded_list)    


    #function
    file_list = glob.glob("unreviewed/*_protein_function_refseq.csv")
    extract_pmid_list(file_list, "evidence", ",", pmid_list, downloaded_list)

    #glycosylation
    file_list = glob.glob("unreviewed/*_proteoform_glycosylation_sites_*.csv")
    extract_pmid_list(file_list, "xref_id", ",", pmid_list, downloaded_list)

    #phosphorylation
    file_list = glob.glob("unreviewed/*_proteoform_phosphorylation_sites_*.csv")
    extract_pmid_list(file_list, "xref_id", ",", pmid_list, downloaded_list)

    #glycosylation
    file_list = glob.glob("unreviewed/*_proteoform_glycation_sites_*.csv")
    extract_pmid_list(file_list, "xref_id", ",", pmid_list, downloaded_list)


    #Reactome
    file_list = glob.glob("unreviewed/*_protein_reactions_reactome.csv")
    extract_pmid_list(file_list, "pmid", ",", pmid_list, downloaded_list)

    #Misc
    file_list = ["downloads/ncbi/medline/misclist.csv"]
    extract_pmid_list(file_list, "pmid", ",", pmid_list, downloaded_list)

    #motif
    file_list = glob.glob("unreviewed/glycan_motif.csv")
    extract_pmid_list(file_list, "pmid", ",", pmid_list, downloaded_list)


    #refseq xref
    extract_pmid_list_from_xref(pmid_list, downloaded_list)

    
    #for pmid in pmid_list:
    #    print pmid
    #sys.exit()

    script_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    n = len(pmid_list)
    n_pass, n_fail = 0, 0
    log_file = "logs/medline-downloads.log"
    with open(log_file, "w") as FL:
        FL.write("Progress\n")

    for pmid in pmid_list:
        time.sleep(1)
        out_file = path_obj["downloads"] + "ncbi/medline/pmid.%s.txt" % (pmid)
        url = "%s?db=pubmed&rettype=medline&id=%s"  % (script_url, pmid)
        res = requests.get(url)
        if res.content.find("PMID") != -1:
            n_pass += 1
            with open(log_file, "a") as FL:
                FL.write("downloaded %s (%s/%s)\n" % (pmid, n_pass, n))
            with open(out_file, 'w') as FW:
                FW.write("%s\n" % (res.content))
        else:
            n_fail += 1
            with open(log_file, "a") as FL:
                FL.write("failed %s (%s/%s)\n" % (pmid, n_fail, n))


    

            
    return


                


if __name__ == '__main__':
        main()
