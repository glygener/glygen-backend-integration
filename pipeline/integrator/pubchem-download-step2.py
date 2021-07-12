import os
import csv
import sys
import json
import commands
import glob
import gzip


sys.path.append('../../glytools/')
import libgly





def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
 
    cid2glytoucan = {}
    data_frame = {}
    for in_file in glob.glob("unreviewed/glycan_xref_pubchem.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[f_list.index("glytoucan_ac")]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            if xref_key == "glycan_xref_pubchem_compound":
                cid2glytoucan[xref_id] = ac

    cmd = "rm -rf " + "downloads/pubchem/compound/current/sdf4glygen/*.sdf"
    x = commands.getoutput(cmd)

    range_dict = {}
    file_list = sorted(glob.glob("downloads/pubchem/compound/sdf/Compound_*.sdf.gz"))
    select_dict = {}
    for in_file in file_list:
        file_name = in_file.split("/")[-1]
        start = int(in_file.split("/")[-1].split(".")[0].split("_")[1])
        end = int(in_file.split("/")[-1].split(".")[0].split("_")[2])
        range_dict[file_name] = {"start":start, "end":end}
    
    for cid in cid2glytoucan:
        cid_int = int(cid)
        for file_name in range_dict:
            cond_list = [cid_int >= range_dict[file_name]["start"]]
            cond_list += [cid_int <= range_dict[file_name]["end"]]
            if False not in cond_list:
                select_dict[file_name] = True
                break

    log_file = "logs/pubchem-downlads-step2.log"
    with open(log_file, "w") as FL:
        FL.write("started logging\n")

    for in_file in file_list:
        file_name = in_file.split("/")[-1]
        if file_name in select_dict:
            cmd = "cp %s %s" % (in_file, "downloads/pubchem/compound/current/sdf4glygen/")
            x = commands.getoutput(cmd)
            with open(log_file, "a") as FL:
                FL.write("Copied %s!\n" % (file_name))



if __name__ == '__main__':
        main()
