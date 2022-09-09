import os,sys
import json
import csv

from optparse import OptionParser

import commands
import glob

import csvutil


__version__="1.0"
__status__ = "Dev"

def get_stats(in_file, rel, seen_dict, c_list, stat_name):

    file_name = in_file.split("/")[-1]
    species = file_name.split("_")[0]

    if species not in seen_dict:
        seen_dict[species] = {}
    if rel not in seen_dict[species]:
        seen_dict[species][rel] = {}

    data_frame = {}
    sep = "," 
    csvutil.load_sheet(data_frame, in_file, [], sep)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        g_type = row[f_list.index("glycosylation_type")].lower()
        saccharide = row[f_list.index("saccharide")].strip()
        if stat_name.find("n_sites") != -1 and g_type.find("n-linked") == -1:
            continue
        if stat_name.find("o_sites") != -1  and g_type.find("o-linked") == -1:
            continue
        if stat_name.find("_with_glycan") != -1 and saccharide == "":
            continue

        val_list = []
        for f in c_list:
            val_list.append(row[f_list.index(f)])
        combo = "^|^".join(val_list)
        seen_dict[species][rel][combo] = True


    return




###############################
def main():

    rel_list = ["1.5.36","1.6.4", "1.7.13","1.8.25","1.9.9","1.10.6","1.11.2", "1.12.3"]

    #rel_list = ["1.11.2", "1.12.3"]
    
    #['uniprotkb_canonical_ac', 'glycosylation_site_uniprotkb', 'amino_acid', 'saccharide', 'glycosylation_type', 'xref_key', 'xref_id', 'src_xref_key', 'src_xref_id']

    c_dict = {
        "glycoproteins":["uniprotkb_canonical_ac"],
        "sites":["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb"],
        "n_sites":["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb"],
        "o_sites":["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb"],
        "sites_with_glycan":["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb"],
        "n_sites_with_glycan":["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb"],
        "o_sites_with_glycan":["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb"]
    }

    source = "all_sources"
    #source = "unicarbkb"

    for stat_name in c_dict:
        seen_dict = {}
        data_dir = "/data/shared/glygen/releases/data/"
        for rel in rel_list:
            reviewed_dir = "/data/shared/glygen/releases/data/v-%s/reviewed/" % (rel)
            pattern = "*_proteoform_glycosylation_sites_*.csv"
            if source != "all_sources":
                pattern = "*_proteoform_glycosylation_sites_*%s*.csv" % (source)
            file_list = glob.glob(reviewed_dir + pattern)
            for in_file in file_list:
                if in_file.find(".stat.csv") != -1:
                    continue
                get_stats(in_file,rel, seen_dict, c_dict[stat_name], stat_name)

        out_file = "logs/stat/%s_%s.csv" % (stat_name, source)
        with open(out_file, "w") as FW:
            row = ["species"]
            for rel in rel_list:
                row.append("release-" + ".".join(rel.split(".")[:-1]))
            FW.write("\"%s\"\n" % ("\",\"".join(row)))
            for species in seen_dict:
                row = [species]
                for rel in rel_list:
                    n = len(seen_dict[species][rel].keys()) if rel in seen_dict[species] else 0
                    row.append(str(n))
                FW.write("\"%s\"\n" % ("\",\"".join(row)))

            


if __name__ == '__main__':
        main()


