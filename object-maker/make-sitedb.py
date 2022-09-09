#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq


import libgly


def add_group_info(order_dict, doc_dict):

    
    win_size = 20
    sorted_dict = sorted(order_dict.items(), key=lambda kv: kv[1])
    site_id_list = []
    for site_id, ordr in sorted_dict:
        site_id_list.append(site_id)

    rel_dict = {}
    for i in range(0, len(site_id_list)):
        site_id_one = site_id_list[i]
        canon, start_pos_one, end_pos_one = site_id_one.split(".")
        cat_list_one = doc_dict[site_id_one]["categories"]
        for j in range(i+1, len(site_id_list)):
            site_id_two = site_id_list[j]
            canon, start_pos_two, end_pos_two = site_id_two.split(".")
            shift = abs(int(end_pos_one) - int(start_pos_two))
            cat_list_two = doc_dict[site_id_two]["categories"]
            if shift <= win_size:
                if site_id_one not in rel_dict:
                    rel_dict[site_id_one] = []
                o_one = {"distance":shift, "direction":"upstream", "id":site_id_two, "categories":cat_list_two}
                rel_dict[site_id_one].append(o_one)
                if site_id_two not in rel_dict:
                    rel_dict[site_id_two] = []
                o_two = {"distance":shift, "direction":"downstream", "id":site_id_one, "categories":cat_list_one}
                rel_dict[site_id_two].append(o_two)

    for site_id in rel_dict:
        doc_dict[site_id]["neighbors"] = rel_dict[site_id]

    return


def main():

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global misc_dir
    global main_dict


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"
    misc_dir = "generated/misc/"


    sec_list = config_obj["sitesections"]

    protein_obj_dict = {}
    record_count = 0

    file_list = glob.glob("jsondb/proteindb/*.json")
    #file_list = glob.glob("jsondb/proteindb/P14210-1.json")
    #file_list = glob.glob("jsondb/proteindb/Q8K4H4*.json")
    #file_list = glob.glob("jsondb/proteindb/Q8R4X3-1*.json")
    #file_list = glob.glob("jsondb/proteindb/Q16555-1*.json")



    for protein_jsonfile in file_list:
        canon = protein_jsonfile.split("/")[-1].split(".")[0]
        doc = json.loads(open(protein_jsonfile,"r").read())
        canon_seq =  doc["sequence"]["sequence"]
        seq_len = len(canon_seq)
        canon = doc["uniprot_canonical_ac"]
        all_sites_dict = {}
        site_dict = {}
        for sec in sec_list:
            if sec in doc:
                for obj in doc[sec]:
                    start_pos, end_pos = -1, -1
                    if "start_pos" in obj:
                        start_pos, end_pos = obj["start_pos"], obj["end_pos"]
                    elif "position" in obj:
                        start_pos, end_pos = obj["position"], obj["position"]
                    if start_pos == -1 or end_pos == -1:
                        continue
                    up_seq, site_seq, down_seq = "", "", ""
                    site_seq = canon_seq[int(start_pos)-1:int(end_pos)]

                    up_seq_start = int(start_pos) - 20 
                    up_seq_start = 1 if up_seq_start < 0 else up_seq_start
                    up_seq_end = int(start_pos) - 1
                    if up_seq_end > 0:
                        up_seq = canon_seq[int(up_seq_start)-1:int(up_seq_end)]
                    down_seq_start = int(end_pos) + 1
                    down_seq_end = int(end_pos) + 20
                    down_seq_start = seq_len if down_seq_start > seq_len else down_seq_start
                    if down_seq_start <= seq_len:
                        down_seq = canon_seq[int(down_seq_start)-1:int(down_seq_end)]

                    site_id = "%s.%s.%s" % (canon,start_pos, end_pos)
                    if site_id not in site_dict:
                        site_dict[site_id] = {
                            "id":site_id, 
                            "start_pos":int(start_pos),
                            "end_pos":int(end_pos),
                            "up_seq":up_seq,
                            "site_seq":site_seq,
                            "down_seq":down_seq,
                            "uniprot_canonical_ac":canon,
                            "species":doc["species"]
                        }

                    o = {}
                    for k in obj:
                        if k not in ["start_pos", "end_pos", "position"]:
                            o[k] = obj[k]

                    if sec not in site_dict[site_id]:
                        site_dict[site_id][sec] = []
                    site_dict[site_id][sec].append(o)
                    xobj = {"type":sec, "start_pos":start_pos, "end_pos":end_pos}
                    x_combo_id = "%s.%s" % (start_pos,end_pos)
                    if sec not in all_sites_dict:
                        all_sites_dict[sec] = {}
                    all_sites_dict[sec][x_combo_id] = True


        all_sites = []
        for site_type in all_sites_dict:
            xobj = {"type": site_type, "site_list": []}
            for r in all_sites_dict[site_type]:
                start_pos, end_pos = r.split(".")
                xobj["site_list"].append({"start_pos": start_pos, "end_pos": end_pos})
            all_sites.append(xobj)
        
       
        doc_dict = {}
        order_dict = {}
        for site_id in site_dict:
            doc = site_dict[site_id]
            doc["all_sites"] = all_sites
            doc["categories"] = []
            empty_sec_count = 0
            for sec in sec_list:
                if sec not in doc:
                    doc[sec] = []
                else:
                    doc["categories"].append(sec + "_flag")
                empty_sec_count += 1 if doc[sec] == [] else 0
            if empty_sec_count == len(sec_list):
                continue
            order_dict[site_id] = int(site_id.split(".")[1])
            doc_dict[site_id] = doc
        add_group_info(order_dict, doc_dict)
       

        for site_id in doc_dict:
            out_file = path_obj["jsondbpath"] + "/sitedb/%s.json" % (site_id)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc_dict[site_id], indent=4)))
            record_count += 1 




    print ("make-sitedb: final created: %s site objects" % (record_count))





if __name__ == '__main__':
    main()

