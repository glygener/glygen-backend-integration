import os,sys
import string
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
from Bio import SeqIO
import csvutil
import commands

__version__="1.0"
__status__ = "Dev"





###############################
def main():


    in_file = "/data/projects/glygen/generated/datasets/compiled/doi_citations.csv"
    data_frame = {}
    csvutil.load_sheet(data_frame, in_file, [], ",")
    f_list = data_frame["fields"]
    seen = {}
    for row in data_frame["data"]:
        doi_id = row[f_list.index("doi_id")]
        s = ""
        for v in row[1:]:
            s += v.strip()
        if s != "":
            seen[doi_id] = True



    reviewed_dir = "/data/projects/glygen/generated/datasets/reviewed/"
    file_list = glob.glob(reviewed_dir + "*_proteoform_glycosylation_sites_*.csv")
    new_dict = {}
    for in_file in file_list:
        data_frame = {}
        csvutil.load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            if xref_key == "protein_xref_doi" and xref_id not in seen:
                new_dict[xref_id] = True


    out_file = "/data/projects/glygen/generated/datasets/compiled/doi_citations4edit.csv"
    with open(out_file, "w") as FW:
        newrow = ["doi_id","title","journal_name","publication_date","authors"]
        FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
        for doi_id in new_dict:
            newrow = [doi_id, "", "", "", ""]
            FW.write("\"%s\"\n" % ("\",\"".join(newrow)))

    cmd = "chmod 775 /data/projects/glygen/generated/datasets/compiled/doi_citations*.csv"
    x = commands.getoutput(cmd)



if __name__ == '__main__':
	main()

