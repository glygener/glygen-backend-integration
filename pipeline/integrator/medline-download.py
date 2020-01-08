import os
import csv
import sys
import json
import commands
import glob
import gzip
import requests


sys.path.append('../../glytools/')
import libgly





def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]
 
    global path_obj
    path_obj = config_obj["pathinfo"]



    pmid_list = {"one":[], "two":[], "three":[], "four":[], "five":[]}
    file_list = glob.glob("unreviewed/*_protein_function_refseq.csv")
    file_list += glob.glob("unreviewed/*_proteoform_glycosylation_sites_unicarbkb.csv")

    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            pmid = row[f_list.index("evidence")]
            if pmid != "":
                pmid_list["one"].append(pmid)
    
    for in_file in [path_obj["downloads"] +  "glytoucan/current/export/pubs.tsv"]:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, "\t")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            pmid = row[f_list.index("PubMedID")]
            if pmid != "":
                pmid_list["two"].append(pmid)


    for in_file in glob.glob("unreviewed/*_protein_reactions_reactome.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            pmid = row[f_list.index("pmid")]
            if pmid != "":
                pmid_list["three"].append(pmid)

    in_file = path_obj["downloads"] + "ncbi/pubmed/medline/misclist.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        pmid = row[f_list.index("pmid")]
        if pmid != "":
            pmid_list["four"].append(pmid)


    #pmid_list = pmid_list["one"] + pmid_list["two"] + pmid_list["three"] + pmid_list["four"]
    
    pmid_list = pmid_list["four"]
    pmid_list = list(set(pmid_list))
    i = 0
    for pmid in pmid_list:
        i += 1
        out_file = path_obj["downloads"] + "ncbi/pubmed/medline/pmid.%s.txt" % (pmid)
        if os.path.isfile(out_file) == False:
            url = "https://www.ncbi.nlm.nih.gov/pubmed/%s?report=medline&format=text" % (pmid)
            res = requests.get(url, verify=False)
            if res.content.strip() != "":
                with open(out_file, 'w') as FW:
                    print "downloaded %s (%s of %s) " % (pmid, i, len(pmid_list))
                    FW.write("%s\n" % (res.content))
        else:
            print "%s already downloaded (%s of %s)" % (pmid, i, len(pmid_list))

    
    return


                


if __name__ == '__main__':
        main()
