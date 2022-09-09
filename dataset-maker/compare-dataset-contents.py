import os,sys
import json
import csv
import time

import commands
import glob

from optparse import OptionParser

import libgly



def load_combo(in_file, selected_fields, seen):
    
    file_name = in_file.split("/")[-1]

    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        if row[f_list.index("xref_key")] not in ["protein_xref_uniprotkb_gly", "protein_xref_pubmed", "protein_xref_doi"]:
            continue
        row[f_list.index("xref_id")] = row[f_list.index("xref_id")].replace("https://doi.org/","")
        
        v_list = [file_name] if "filename" in selected_fields else [] 
        for f in selected_fields:
            v_list.append(row[f_list.index(f)] if f in f_list else "N/A")
        combo = "^|^".join(v_list)
        seen[combo] = True

    return



def main():


    file_list = glob.glob("unreviewed/*_proteoform_glycosylation_sites_*.csv")
    #selected_fields = ["amino_acid", "start_pos", "end_pos", "start_aa", "end_aa"]
    selected_fields = ["amino_acid", "start_pos", "end_pos","filename"]
    for in_file in file_list:
        seen = {}
        load_combo(in_file, selected_fields, seen)
        for combo in seen:
            parts = combo.split("^|^")
            #if parts[2] == parts[3] and parts[1].strip() == "":
            if parts[2] == "" or parts[3] == "":
                print combo



    exit()


   
    old_ver = "1.11.2"
    file_list = glob.glob("unreviewed/*_proteoform_glycosylation_sites_*.csv")
    rel_dir = "/data/shared/glygen/releases/data/v-1.11.2/"
    selected_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","glycosylation_type",
            "saccharide","xref_key","xref_id"]

    log_dict = {"new":{}, "missing":{}}
    for in_file_two in file_list:
        file_name = in_file_two.split("/")[-1]
        in_file_one = rel_dir + "/reviewed/" + file_name
        if os.path.isfile(in_file_one) == False:
            continue
        seen = {"one":{}, "two":{}}
        load_combo(in_file_one, selected_fields, seen["one"])
        load_combo(in_file_two, selected_fields, seen["two"])
        for combo in seen["one"]:
            if combo not in seen["two"]:
                if file_name not in log_dict["missing"]:
                    log_dict["missing"][file_name] = []
                log_dict["missing"][file_name].append(combo.split("^|^"))
        
        for combo in seen["two"]:
            if combo not in seen["one"]:
                if file_name not in log_dict["new"]:
                    log_dict["new"][file_name] = []
                log_dict["new"][file_name].append(combo.split("^|^"))

    cmd = "rm -f logs/misc/missing_in_current/*"
    x = commands.getoutput(cmd)
    cmd = "rm -f logs/misc/new_in_current/*"
    x = commands.getoutput(cmd)

    out_dir = "logs/misc/new_in_current/"
    hrow = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","saccharide"]
    hrow += ["xref_key","xref_id"]
    for file_name in log_dict["new"]:
        with open(out_dir + "%s" % (file_name), "w") as FW:
            FW.write("\"%s\"\n" % ("\",\"".join(hrow)))
            for row in log_dict["new"][file_name]:
                FW.write("\"%s\"\n" % ("\",\"".join(row)))
    
    out_dir = "logs/misc/missing_in_current/"
    hrow = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","saccharide"]
    hrow += ["xref_key","xref_id"]
    for file_name in log_dict["missing"]:
        with open(out_dir + "%s" % (file_name), "w") as FW:
            FW.write("\"%s\"\n" % ("\",\"".join(hrow)))
            for row in log_dict["missing"][file_name]:
                FW.write("\"%s\"\n" % ("\",\"".join(row)))
 
           
 
if __name__ == '__main__':
            
    main()



