import os
import sys
import json
import glob
from optparse import OptionParser

import libgly


def load_is_conjugated():
   
    file_list = glob.glob("outdir/*_glycan2protein.csv")
    for in_file in file_list:
        data_frame = {}
        delim = "," if in_file.split(".")[-1] == "csv" else "\t"
        libgly.load_sheet(data_frame, in_file, delim)
        for row in data_frame["data"]:
            is_conjugated[row[0]] = True
    
    return

def load_glycan2composition(in_file):

    seen = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[0]
        if glytoucan_ac not in seen:
            seen[glytoucan_ac] = {}
        for f in f_list[1:-1]:
            seen[glytoucan_ac][f] = str(int(row[f_list.index(f)].replace("+", "")))

    return seen

def load_glycan2protein(species, gtype_list):

    file_list = glob.glob("unreviewed/%s_proteoform_glycosylation_sites_*.csv"%(species))
    seen = {}
    for in_file in file_list:
        data_frame = {}
        delim = "," if in_file.split(".")[-1] == "csv" else "\t"
        libgly.load_sheet(data_frame, in_file, delim)
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            glytoucan_ac = row[f_list.index("saccharide")]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
            g_type = row[f_list.index("glycosylation_type")]
            if g_type not in gtype_list:
                continue
            if xref_key != "protein_xref_pubmed":
                continue
            if canon == "" or glytoucan_ac == "":
                continue
            if glytoucan_ac not in seen:
                seen[glytoucan_ac] = {}
            seen[glytoucan_ac][canon] = True

    return seen

def load_simple_matrix(in_file,row_field, col_field):

    tmp_dict = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        row_id = row[f_list.index(row_field)].strip()
        col_id = row[f_list.index(col_field)].strip()
        if row_id == "" or col_id == "":
            continue
        if row_field == "glytoucan_ac" and row_id not in is_conjugated:
            continue
        if col_field == "glytoucan_ac" and col_id not in is_conjugated:
            continue
        if row_id not in tmp_dict:
            tmp_dict[row_id] = {}
        tmp_dict[row_id][col_id] = True

    return tmp_dict


def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " )
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")
    parser.add_option("-d","--dataset",action="store",dest="dataset",help="glycan2protein")

    (options,args) = parser.parse_args()
    for file in ([options.species, options.dataset]):
        if not (file):
            parser.print_help()
            sys.exit(0)
    
    
    global is_conjugated



    species = options.species
    dataset = options.dataset

    is_conjugated = {}
    load_is_conjugated()
   

    data_matrix = {}
    out_file = ""
    if dataset == "glycan2protein":
        data_matrix = load_glycan2protein(species, ["N-linked"])
        out_file = "outdir/%s_%s.csv" % (species, dataset)
    elif dataset == "glycan2composition":
        in_file = "unreviewed/glycan_monosaccharide_composition.csv"
        data_matrix = load_glycan2composition(in_file)
        out_file = "outdir/%s.csv" % (dataset)
    elif dataset == "glycan2compositionadv":
        in_file = "unreviewed/glycan_monosaccharide_composition_advanced.csv"
        data_matrix = load_glycan2composition(in_file)
        out_file = "outdir/%s.csv" % (dataset)
    elif dataset == "motif2glycan":
        in_file = "unreviewed/glycan_motif.csv"
        row_field, col_field = "motif_ac", "glytoucan_ac"
        data_matrix = load_simple_matrix(in_file,row_field, col_field)
        out_file = "outdir/%s.csv" % (dataset)
    elif dataset == "enzyme2glycan":
        in_file = "unreviewed/glycan_enzyme.csv"
        row_field, col_field = "uniprotkb_canonical_ac", "glytoucan_ac"
        data_matrix = load_simple_matrix(in_file,row_field, col_field)
        out_file = "outdir/%s.csv" % (dataset)


    FW = open(out_file, "w")
    row = ["row_id", "col_id", "value"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))
    for row_id in data_matrix:
        for col_id in data_matrix[row_id]:
            row = [row_id, col_id, data_matrix[row_id][col_id]]
            FW.write("\"%s\"\n" % ("\",\"".join(row)))
    FW.close()
            
    return


                


if __name__ == '__main__':
        main()
