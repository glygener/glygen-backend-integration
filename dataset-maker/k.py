import os
import csv
import sys
import json
import glob
import subprocess
import libgly

def main():

    ac2canon = {}
    file_list = glob.glob("unreviewed/*_protein_masterlist.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            ac = canon.split("-")[0]
            ac2canon[ac] = canon
            for f in ["reviewed_isoforms","unreviewed_isoforms"]:
                ac = row[f_list.index(f)].split("-")[0]
                ac2canon[ac] = canon



    omaid2ac = {}
    in_file = "downloads/oma/current/oma-uniprot.txt"
    with open(in_file, "r") as FR:
        for line in FR:
            if line[0] == "#":
                continue
            row = line.strip().split("\t")
            oma_id, ac = row[0], row[1]
            if ac not in ac2canon:
                continue
            if oma_id not in omaid2ac:
                omaid2ac[oma_id] = {}
            omaid2ac[oma_id][ac] = True

    cls_dict = {}
    file_list = glob.glob("downloads/oma/current/*.json")
    for in_file in file_list:
        doc = json.loads(open(in_file, "r").read())
        for obj in doc:
            o_one, o_two = obj["entry_1"], obj["entry_2"]
            oma_id_one, canonical_id_one = o_one["omaid"], o_one["canonicalid"]
            tax_id_one = str(o_one["species"]["taxon_id"])
            oma_group_one = str(o_one["oma_group"])
            oma_id_two, canonical_id_two = o_two["omaid"], o_two["canonicalid"]
            tax_id_two = str(o_two["species"]["taxon_id"])
            oma_group_two = str(o_two["oma_group"])
            ac_list_one = omaid2ac[oma_id_one].keys() if oma_id_one in omaid2ac else []
            for ac in ac_list_one:
                if ac not in ac2canon:
                    continue
                canon = ac2canon[ac]
                if oma_group_one not in cls_dict:
                    cls_dict[oma_group_one] = {}
                cls_dict[oma_group_one][canon] = True
            ac_list_two = omaid2ac[oma_id_two].keys() if oma_id_two in omaid2ac else []
            for ac in ac_list_two:
                if ac not in ac2canon:
                    continue
                canon = ac2canon[ac]
                if oma_group_two not in cls_dict:
                    cls_dict[oma_group_two] = {}
                cls_dict[oma_group_two][canon] = True

    for cls in cls_dict:
        n = len(cls_dict[cls].keys())
        print n, cls

if __name__ == '__main__':
        main()
