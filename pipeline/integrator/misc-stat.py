import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient

sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"



###############################
def main():


    site_dict = {"glcoyslation":{}, "mutation":{}}
    site_type_dict = {"glcoyslation":{}, "mutation":{} }
    for in_file in glob.glob("unreviewed/*_protein_glycosylation_motifs.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            start_pos = row[f_list.index("start_pos")]
            motif = row[f_list.index("motif")]
            motif_type = motif[0]
            if canon not in site_dict["glcoyslation"]:
                site_dict["glcoyslation"][canon] = []
                site_type_dict["glcoyslation"][canon] = []
            if start_pos not in site_dict["glcoyslation"][canon]:
                site_dict["glcoyslation"][canon].append(start_pos)
                site_type_dict["glcoyslation"][canon].append(motif_type)

    for in_file in glob.glob("unreviewed/*_protein_mutation.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            start_pos = row[f_list.index("aa_pos")]
            if canon not in site_dict["mutation"]:
                site_dict["mutation"][canon] = []
            if start_pos not in site_dict["mutation"][canon]:
                site_dict["mutation"][canon].append(start_pos)

    newrow = ["uniprotkb_canonical_ac","glycosylation_site_count", "mutation_site_count","motif_types"]
    print "\"%s\""  % ("\",\"".join(newrow))

    canon_list = sorted(list(set(site_dict["glcoyslation"].keys() + site_dict["mutation"].keys())))
    for canon in canon_list:
        newrow = [canon]
        g1, g2, m1 = "", "", ""
        if canon in site_dict["glcoyslation"]:
            g1 = str(len(site_dict["glcoyslation"][canon]))
            g2 = ";".join(sorted(list(set(site_type_dict["glcoyslation"][canon]))))
        if canon in site_dict["mutation"]:
            m1 = str(len(site_dict["mutation"][canon]))
        newrow = [canon, g1, m1, g2]
        print "\"%s\""  % ("\",\"".join(newrow))


if __name__ == '__main__':
	main()

