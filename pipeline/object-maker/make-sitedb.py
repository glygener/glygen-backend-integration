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
import commands


sys.path.append('../../glytools/')
import libgly


def main():

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global misc_dir
    global main_dict


    config_file = "../../conf/config-1.1.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"
    misc_dir = "generated/misc/"


    sec_list = [
        "glycosylation"
        ,"mutagenesis"
        ,"snv"
        ,"site_annotation"
        ,"phosphorylation"
        ,"glycation"
    ]
    protein_obj_dict = {}
    record_count = 0
    file_list = glob.glob("jsondb/proteindb/*.json")
    #file_list = glob.glob("jsondb/proteindb/P1421*.json")
    #file_list = glob.glob("jsondb/proteindb/Q8K4H4*.json")
    #file_list = glob.glob("jsondb/proteindb/Q8R4X3-1*.json")


    for protein_jsonfile in file_list:
        canon = protein_jsonfile.split("/")[-1].split(".")[0]
        doc = json.loads(open(protein_jsonfile,"r").read())
        canon_seq =  doc["sequence"]["sequence"]
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
                    site_seq = canon_seq[int(start_pos)-1:int(end_pos)]
                    site_id = "%s.%s.%s" % (canon,start_pos, end_pos)
                    if site_id not in site_dict:
                        site_dict[site_id] = {
                            "id":site_id, 
                            "start_pos":int(start_pos),
                            "end_pos":int(end_pos),
                            "site_seq":site_seq,
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


        # Ingest single position section data of a site
        # to parent range site section 
        for site_id in site_dict:
            canon,start_pos, end_pos = site_id.split(".")
            if start_pos != end_pos:
                for i in xrange(int(start_pos), int(end_pos) + 1):
                    in_site_id = "%s.%s.%s" % (canon, i, i)
                    if in_site_id in site_dict:
                        for sec in sec_list:
                            if sec not in site_dict[in_site_id]:
                                continue
                            if sec not in site_dict[site_id]:
                                site_dict[site_id][sec] = []
                            site_dict[site_id][sec] += site_dict[in_site_id][sec]

        all_sites = []
        for site_type in all_sites_dict:
            xobj = {"type": site_type, "site_list": []}
            for r in all_sites_dict[site_type]:
                start_pos, end_pos = r.split(".")
                xobj["site_list"].append({"start_pos": start_pos, "end_pos": end_pos})
            all_sites.append(xobj)

        for site_id in site_dict:
            doc = site_dict[site_id]
            doc["all_sites"] = all_sites
            empty_sec_count = 0
            for sec in sec_list:
                if sec not in doc:
                    doc[sec] = []
                empty_sec_count += 1 if doc[sec] == [] else 0
            if empty_sec_count == len(sec_list):
                #print site_id, empty_sec_count, len(sec_list)
                continue

            out_file = path_obj["jsondbpath"] + "/sitedb/%s.json" % (site_id)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            record_count += 1 
            #if record_count%1000 == 0:
            #    print " ... created %s objects" % (record_count)
    
    print " ... final created: %s site objects" % (record_count)




if __name__ == '__main__':
    main()

