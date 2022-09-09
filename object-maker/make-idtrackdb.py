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

import datetime
import pytz

import libgly

import subprocess





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

    data_release_dir = "/data/shared/glygen/releases/data/"
    out_json = {}
    
    url = "https://api.glygen.org//misc/verlist/"
    cmd = "curl -s -k %s" % (url)
    res = subprocess.getoutput(cmd)
    release_list = sorted(json.loads(res)) + ["current"]


    species_obj = {}
    in_file = path_obj["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)

    ref_tax_id_list = []
    for k in species_obj:
        if species_obj[k]["is_reference"] == "yes":
            if str(species_obj[k]["tax_id"]) not in ref_tax_id_list:
                ref_tax_id_list.append(str(species_obj[k]["tax_id"]))

    ts_dict = {}
    pattern = "/data/shared/glygen/releases/data/v-*/reviewed/release-notes.txt"
    for release_file in glob.glob(pattern):
        rel = release_file.split("/")[-3].split("-")[-1]
        if rel not in release_list:
            continue

        mm,dd,yy = open(release_file, "r").read().split("\n")[0].split(" ")[1].split("/")
        ts_dict[rel] = "%s-%s-%s 00:00:00" % (yy, mm, dd)
    ts = datetime.datetime.now(pytz.timezone('US/Eastern')).strftime('%Y-%m-%d')
    ts_dict["current"] = ts + " 00:00:00"


    history_dict = {}
    in_file = "reviewed/protein_uniprotkb_accession_history.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        old_ac = row[f_list.index("uniprotkb_ac_old")]
        current_ac = row[f_list.index("uniprotkb_ac_current")]
        if old_ac not in history_dict:
            history_dict[old_ac] = []
        if current_ac not in history_dict[old_ac]:
            history_dict[old_ac].append(current_ac)


    in_file = "reviewed/glycan_glytoucan_accession_history.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        old_ac = row[f_list.index("glytoucan_ac_old")]
        current_ac = row[f_list.index("glytoucan_ac_current")]
        if old_ac not in history_dict:
            history_dict[old_ac] = []
        if current_ac not in history_dict[old_ac]:
            history_dict[old_ac].append(current_ac)




    ac2canon = {}
    record_type_dict = {}
    track_dict = {}
    seen_taxid = {}
    for rel in release_list:
        reviewed_dir = data_release_dir + "v-%s/reviewed/" % (rel)
        if rel == "current":
            reviewed_dir = "reviewed/"

        seen_taxid[rel] = {}
        protein_file_list, glycan_file_list, glycan_species_file_list = [], [], []
        if rel in ["1.0.13", "1.4.5"]:
            pattern = reviewed_dir + "*_protein_idmapping.csv" 
            protein_file_list += glob.glob(pattern)
            pattern = reviewed_dir + "*_glycan_idmapping.csv" 
            glycan_file_list += glob.glob(pattern)
            if glycan_file_list == []:
                pattern = reviewed_dir + "*_glycan_properties.csv" 
                glycan_file_list += glob.glob(pattern)
            glycan_species_file_list = glycan_file_list
        else:
            pattern = reviewed_dir + "*_protein_masterlist.csv"
            protein_file_list += glob.glob(pattern)
            pattern = reviewed_dir + "glycan_masterlist.csv" 
            glycan_file_list += glob.glob(pattern)
            pattern = reviewed_dir + "glycan_taxonomy.csv" 
            glycan_species_file_list += glob.glob(pattern)
            pattern = reviewed_dir + "glycan_species.csv" 
            glycan_species_file_list += glob.glob(pattern)

        for in_file in protein_file_list:
            mol = in_file.split("/")[-1].split("_")[1]
            species = in_file.split("/")[-1].split("_")[0]
            tax_id = str(species_obj[species]["tax_id"])
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            for row in data_frame["data"]:
                main_id = row[0]
                if main_id not in seen_taxid:
                    seen_taxid[main_id] = {}
                seen_taxid[main_id][tax_id] = True
                if main_id not in track_dict:
                    track_dict[main_id] = {}
                track_dict[main_id][rel] = True
                record_type_dict[main_id] = mol
                if rel not in ac2canon:
                    ac2canon[rel] = {}
                ac = main_id.split("-")[0]
                ac2canon[rel][ac] = main_id
        for in_file in glycan_file_list:
            mol = in_file.split("/")[-1].split("_")[0]
            if rel in ["1.0.13", "1.4.5"]:
                mol = in_file.split("/")[-1].split("_")[1]
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            for row in data_frame["data"]:
                main_id = row[0]
                if main_id not in track_dict:
                    track_dict[main_id] = {}
                track_dict[main_id][rel] = True
                record_type_dict[main_id] = mol
       
 
        for in_file in glycan_species_file_list:
            tax_id = ""
            if rel in ["1.0.13", "1.4.5"]:
                species = in_file.split("/")[-1].split("_")[0]
                tax_id = str(species_obj[species]["tax_id"])
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                main_id = row[0]
                if rel not in ["1.0.13", "1.4.5"]:
                    tax_id = row[f_list.index("tax_id")]
                if tax_id not in species_obj:
                    continue
                if species_obj[tax_id]["is_reference"] != "yes":
                    continue
                if main_id not in seen_taxid[rel]:
                    seen_taxid[rel][main_id] = {}
                seen_taxid[rel][main_id][tax_id] = True    
    


    record_stat = {}
    record_count = 0
    msg_dict = {}
    for main_id in track_dict:
        record_type = record_type_dict[main_id]
        prev_status = False
        if record_type not in msg_dict:
            msg_dict[record_type]  = {}
        for rel in release_list:
            status = rel in track_dict[main_id]
            replacement_id_list = []
            replacement_status = False
            if record_type == "protein":
                old_ac = main_id.split("-")[0]
                if old_ac in history_dict:
                    for new_ac in history_dict[old_ac]:
                        if new_ac in ac2canon[rel]:
                            replacement_id = ac2canon[rel][new_ac]
                            replacement_status = "yes" if rel in track_dict[replacement_id] else "no"
                            combo_value = "%s %s" % (replacement_id, replacement_status)
                            replacement_id_list.append(combo_value)

            tax_id_list = ["unknown"]
            if record_type == "protein":
                if main_id in seen_taxid:
                    tax_id_list = seen_taxid[main_id].keys()
            else:
                if main_id in seen_taxid[rel]:
                    tax_id_list = seen_taxid[rel][main_id].keys()
            o_list = []
            if status == False and prev_status == True:
                msg = "discontinued in data release %s" % (rel)
                o_list.append({"type":"removed", "timestamp":ts_dict[rel],"description":msg})
                if replacement_id_list !=  []:
                    for combo_value in replacement_id_list:
                        replacement_id, replacement_status = combo_value.split(" ")
                        if replacement_status == "yes":
                            msg = "replaced by %s in data release %s" % (replacement_id, rel)
                            o_list.append({"type":"replaced", "replacement_id":replacement_id,
                                "timestamp":ts_dict[rel],"description":msg})

                if rel not in record_stat:
                    record_stat[rel] = {}
                if record_type not in record_stat[rel]:
                    record_stat[rel][record_type] = {}
                for tax_id in tax_id_list:
                    if tax_id not in record_stat[rel][record_type]:
                        record_stat[rel][record_type][tax_id] = {"removed":0, "added":0}
                    record_stat[rel][record_type][tax_id]["removed"] += 1

            if status == True and prev_status == False:
                msg = "introduced in data release %s" % (rel)
                o_list.append({"type":"added", "timestamp":ts_dict[rel],"description":msg})
                if rel not in record_stat:
                    record_stat[rel] = {}
                if record_type not in record_stat[rel]:
                    record_stat[rel][record_type] = {}
                for tax_id in tax_id_list:
                    if tax_id not in record_stat[rel][record_type]:
                        record_stat[rel][record_type][tax_id] = {"removed":0, "added":0}
                    record_stat[rel][record_type][tax_id]["added"] += 1

            prev_status = status
            if o_list != []:
                if main_id not in msg_dict[record_type]:
                    msg_dict[record_type][main_id] = []
                for o in o_list:
                    msg_dict[record_type][main_id].append(o)
                    if "replacement_id" in o:
                        replacor_id = o["replacement_id"]
                        msg = "inherits %s as of data release %s" % (main_id, rel)
                        oo = {"type":"inherits", "inherited_id":main_id,
                                "timestamp":ts_dict[rel],"description":msg}
                        if replacor_id not in msg_dict[record_type]:
                            msg_dict[record_type][replacor_id] = []
                        msg_dict[record_type][replacor_id].append(oo)




    for record_type in msg_dict:
        for main_id in msg_dict[record_type]:
            doc = {"record_id":main_id,"history":msg_dict[record_type][main_id]}
            out_file = path_obj["jsondbpath"] + "/idtrackdb/%s.json" % (main_id)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            record_count += 1 


    out_json = {}
    for rel in release_list:
        for record_type in ["protein", "glycan"]:
            for tax_id in sorted(ref_tax_id_list) + ["unknown"]:
                tax_name = species_obj[tax_id]["long_name"] if tax_id in species_obj else tax_id
                added, removed = 0, 0
                if rel in record_stat:
                    if record_type in record_stat[rel]:
                        if tax_id in record_stat[rel][record_type]:
                            added = record_stat[rel][record_type][tax_id]["added"]
                            removed = record_stat[rel][record_type][tax_id]["removed"]
                o = {"release":rel, "release_date":ts_dict[rel], "added":added, "removed":removed}
                if tax_name not in out_json:
                    out_json[tax_name] = {}
                if record_type not in out_json[tax_name]:
                    out_json[tax_name][record_type] = []
                out_json[tax_name][record_type].append(o)

    out_file = path_obj["jsondbpath"] + "/logs/idtrackdb.json" 
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(out_json, indent=4)))


    print ("make-idtrackdb: final created: %s idtrack objects" % (record_count))





if __name__ == '__main__':
    main()

