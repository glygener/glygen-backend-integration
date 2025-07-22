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
import time
import libgly
import csvutil

import subprocess
import requests


def output_glygen(msg_dict):

    record_count = 0
    for record_type in msg_dict:
        for main_id in msg_dict[record_type]:
            doc = {"record_id":main_id,"recordtype":record_type,
                    "history":msg_dict[record_type][main_id]}
            out_file =  "jsondb/idtrackdb/%s.glygen.%s.json" % (record_type,main_id)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            record_count += 1 

    return record_count

def output_batches(batch_doc_list):

    record_count = 0
    for doc in batch_doc_list:
        idx = doc["batchidx"]
        out_file = "jsondb/idtrackdb/protein.batch.%s.json" % (idx)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        record_count += 1

    return record_count



def output_glytoucan_only(history_dict, glygen_ac_dict, glytoucan_ac_dict):

    record_count = 0
    #accessions in glytoucan but not in glygen
    for ac in glytoucan_ac_dict:
        if ac not in glygen_ac_dict:
            out_file = "jsondb/idtrackdb/glycan.glytoucan.%s.json" % (ac)
            doc = {"record_id":ac, "history":[], "recordtype":"glycan"}
            o = {}
            o["type"] = "never_in_glygen_current_in_glytoucan"
            o["description"] = "GlyTouCan accession never been in GlyGen"
            if glytoucan_ac_dict[ac].find("Discontinued") != -1:
                o["type"] = "never_in_glygen_discontinued_in_glytoucan"
                o["description"] = "Discontinued GlyTouCan accession never been in GlyGen" 
            if ac in history_dict:
                new_ac = history_dict[ac][0]
                o["type"], o["replacement_id_list"] = "replacement_not_in_glygen", [new_ac]
                o["description"] = "Old GlyToucan accession replaced by %s" % (new_ac)
                if new_ac in glygen_ac_dict:
                    o["type"], o["replacement_id_list"] = "replacement_in_glygen", [new_ac]
            doc["history"].append(o)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            record_count += 1

    return record_count


def output_uniprotkb_only(history_dict, glygen_ac_dict, ac2canon):

    record_count = 0
    ac_list = list(history_dict.keys())
    for ac in ac_list:
        if ac not in glygen_ac_dict:
            if ac in history_dict:
                canon_list = []
                for new_ac in history_dict[ac]:
                    if new_ac in glygen_ac_dict and new_ac in ac2canon["current"]:
                        canon_list.append(ac2canon["current"][new_ac])
                if canon_list != []:
                    r = ", ".join(canon_list)
                    o = {
                        "description":"Old UniProtKB accession replaced by %s" % (r),
                        "type":"replacement_in_glygen",
                        "replacement_id_list":[canon_list]
                    }
                    doc = {"record_id":ac, "history":[o], "recordtype":"protein"}
                    out_file = "jsondb/idtrackdb/protein.uniprotkb.%s.json" % (ac)
                    with open(out_file, "w") as FW:
                        FW.write("%s\n" % (json.dumps(doc, indent=4)))
                    record_count += 1

    return record_count



def sort_release_list(tmp_list, reversed_flag):

    factor_list = [100000000, 1000, 1]
    rel_dict = {}
    for rel in tmp_list:
        parts = rel.split(".")
        ordr = 0
        for i in range(0,len(parts)):
            ordr += factor_list[i]*int(parts[i])
        rel_dict[ordr] = rel
    
    release_list = []

    for ordr in sorted(rel_dict, reverse=reversed_flag):
        release_list.append(rel_dict[ordr])


    seen = {}
    filtered_release_list = []
    for r in release_list:
        rr = r.split(".")[0] + "." + r.split(".")[1]
        if rr not in seen:
            seen[rr] = True
            filtered_release_list.append(r)


    return filtered_release_list








def main():

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global main_dict


    #DEBUG = True
    DEBUG = False

    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"

    data_release_dir = "/data/shared/glygen/releases/data/"
    out_json = {}
    
    url = "https://api.glygen.org/misc/verlist"
    res = requests.get(url,  verify=False)
    res_obj = json.loads(res.content)
    release_list = sort_release_list(res_obj, False)
    release_list = release_list + ["current"]

    #print (release_list)
    #exit()

    ts_dict = {}
    file_list = glob.glob("/data/shared/glygen/releases/data/v-*/reviewed/human_protein_masterlist.csv")
    file_list += glob.glob("/data/shared/glygen/releases/data/v-*/reviewed/human_protein_idmapping.csv")
    for in_file in file_list:
        rel = in_file.split("/")[6].split("-")[-1]
        if rel not in release_list:
            continue
        ts = datetime.datetime.fromtimestamp(os.path.getmtime(in_file)).strftime('%Y-%m-%d')
        ts_dict[rel] = ts + " 00:00:00"
    ts = datetime.datetime.now(pytz.timezone('US/Eastern')).strftime('%Y-%m-%d')
    ts_dict["current"] = ts + " 00:00:00"




    all_ac_file = "downloads/ebi/current/all_accessions.csv"
    if DEBUG:
        all_ac_file = "downloads/ebi/current/toy.csv"
        release_list = release_list[-3:]


    species_obj = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)

    ref_tax_id_list = []
    for k in species_obj:
        if species_obj[k]["is_reference"] == "yes":
            if str(species_obj[k]["tax_id"]) not in ref_tax_id_list:
                ref_tax_id_list.append(str(species_obj[k]["tax_id"]))




    log_file = "logs/make-idtrackdb.log"
    msg = "make-idtrackdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")


    msg = "make-idtrackdb: loading all_accessions.csv"
    csvutil.write_log_msg(log_file, msg, "a")

    history_dict = {}
    bucket_dict = {"primary":"", "secondary":"", "discontinued":""}
    batch_idx = 0
    batch_doc_list = []
    buffer_size = 100000
    with open(all_ac_file, "r") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            if lcount == 1:
                continue
            #ac, status = line.strip().replace(" ", "").split(",")
            ac, status, primary_ac = line.strip().replace(" ", "").split(",")
            if status == "primary":
                bucket_dict["primary"] += "%s " % (ac)
                if len(bucket_dict["primary"]) > buffer_size:
                    batch_idx += 1
                    aclist = bucket_dict["primary"].strip().replace(" ", ",")
                    doc = {"batchidx":batch_idx, "accessions":aclist,"status":"primary"}
                    doc["recordtype"] = "protein"
                    batch_doc_list.append(doc)
                    bucket_dict["primary"] = ""
            elif status == "secondary" and primary_ac != "":
                bucket_dict["secondary"] += "%s:%s " % (ac, primary_ac)
                if len(bucket_dict["secondary"]) > buffer_size:
                    batch_idx += 1
                    aclist = bucket_dict["secondary"].strip().replace(" ", ",")
                    doc = {"batchidx":batch_idx, "accessions":aclist,"status":"secondary"}
                    doc["recordtype"] = "protein"
                    batch_doc_list.append(doc)
                    bucket_dict["secondary"] = ""
                if ac not in history_dict:
                    history_dict[ac] = []
                if primary_ac not in history_dict[ac]:
                    history_dict[ac].append(primary_ac)
            elif status == "obsolete":
                bucket_dict["discontinued"] += "%s " % (ac)
                if len(bucket_dict["discontinued"]) > buffer_size:
                    batch_idx += 1
                    aclist = bucket_dict["discontinued"].strip().replace(" ", ",")
                    doc = {"batchidx":batch_idx, "accessions":aclist, "status":"discontinued"}
                    doc["recordtype"] = "protein"
                    batch_doc_list.append(doc)
                    bucket_dict["discontinued"] = ""

    msg = "make-idtrackdb: loading glytoucan_ac_dict"
    csvutil.write_log_msg(log_file, msg, "a")

    glytoucan_ac_dict = {}
    in_file = "reviewed/glycan_glytoucanidlist.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac  = row[f_list.index("glytoucan_ac")]
        if ac.strip() != "":
            glytoucan_ac_dict[ac] = "GlyTouCan accession"
            status, new_ac = row[f_list.index("status")], row[f_list.index("replacement_ac")]
            if status == "discontinued":
                glytoucan_ac_dict[ac] = "Discontinued GlyTouCan accession"
            elif status == "replaced":
                glytoucan_ac_dict[ac] = "Old GlyTouCan accession "
                glytoucan_ac_dict[ac] += "replaced by %s" % (new_ac)

    in_file = "reviewed/glycan_glytoucan_accession_history.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        old_ac = row[f_list.index("glytoucan_ac_old")]
        current_ac = row[f_list.index("glytoucan_ac_current")]
        if current_ac.strip() != "":
            glytoucan_ac_dict[current_ac] = "GlyTouCan accession"
            if old_ac not in history_dict:
                history_dict[old_ac] = []
            if current_ac not in history_dict[old_ac]:
                history_dict[old_ac].append(current_ac)
            glytoucan_ac_dict[old_ac] = "Old GlyTouCan accession "
            glytoucan_ac_dict[old_ac] += "replaced by %s" % (current_ac)
        else:
            glytoucan_ac_dict[old_ac] = "Discontinued GlyTouCan accession"



    msg = "make-idtrackdb: loading track_dict"
    csvutil.write_log_msg(log_file, msg, "a")
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
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                main_id = row[f_list.index("uniprotkb_canonical_ac")]
                isoform_list = [
                    row[f_list.index("reviewed_isoforms")],
                    row[f_list.index("unreviewed_isoforms")]
                ]

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
                for isoform in isoform_list:
                    isoform_ac = isoform.split("-")[0]
                    ac2canon[rel][isoform_ac] = main_id

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

    main_id_list = list(track_dict.keys())

    record_stat = {}
    msg_dict = {}
    glygen_ac_dict = {}
    for main_id in main_id_list:
        glygen_ac_dict[main_id.split("-")[0]] = True
        record_type = record_type_dict[main_id]
        prev_status = False
        if record_type not in msg_dict:
            msg_dict[record_type]  = {}
        for rel in release_list:
            status = rel in track_dict[main_id]
            rep_id_list = []
            replacement_status = False
            if record_type in ["protein"]:
                old_ac = main_id.split("-")[0]
                canon_list = []
                if old_ac in history_dict:
                    for new_ac in history_dict[old_ac]:
                        if new_ac in ac2canon[rel]:
                            canon = ac2canon[rel][new_ac]
                            if rel in track_dict[canon]:
                                rep_id_list.append(canon)
            elif record_type in ["glycan"]:
                if main_id in history_dict:
                    for gtc_ac in history_dict[main_id]:
                        replacement_status = "no" 
                        if gtc_ac in track_dict:
                            if rel in track_dict[gtc_ac]:
                                rep_id_list.append(gtc_ac)
            tax_id_list = ["unknown"]
            if record_type == "protein":
                if main_id in seen_taxid:
                    tax_id_list = seen_taxid[main_id].keys()
            else:
                if main_id in seen_taxid[rel]:
                    tax_id_list = seen_taxid[rel][main_id].keys()
            o_list = []
            xref_status = ""
            xref_status = "Old " if main_id in history_dict else xref_status
            xref_status = "Old " if main_id.split("-")[0] in history_dict else xref_status
            if status == False and prev_status == True:
                msg_str = "%s%s accession replaced by %s"
                resource = "UniProtKB" if record_type == "protein" else "GlyTouCan"
                if rep_id_list !=  []:
                    o_list.append({
                        "type":"replacement_in_glygen", 
                        "replacement_id_list":rep_id_list,
                        "description":msg_str % (xref_status,resource,",".join(rep_id_list))
                    })
                if o_list == []:
                    msg = "%s%s accession discontinued in GlyGen in release %s" % (xref_status,resource, rel)
                    o_list.append({"type":"discontinued_in_glygen", "description":msg})
                if rel not in record_stat:
                    record_stat[rel] = {}
                if record_type not in record_stat[rel]:
                    record_stat[rel][record_type] = {}
                for tax_id in tax_id_list:
                    if tax_id not in record_stat[rel][record_type]:
                        record_stat[rel][record_type][tax_id] = {"discontinued_in_glygen":0, "added":0}
                    record_stat[rel][record_type][tax_id]["discontinued_in_glygen"] += 1

            if status == True and prev_status == False:
                msg = "introduced in data release %s" % (rel)
                o_list.append({"type":"added", "description":msg})
                if rel not in record_stat:
                    record_stat[rel] = {}
                if record_type not in record_stat[rel]:
                    record_stat[rel][record_type] = {}
                for tax_id in tax_id_list:
                    if tax_id not in record_stat[rel][record_type]:
                        record_stat[rel][record_type][tax_id] = {"discontinued_in_glygen":0, "added":0}
                    record_stat[rel][record_type][tax_id]["added"] += 1

            prev_status = status
            if o_list != []:
                if main_id not in msg_dict[record_type]:
                    msg_dict[record_type][main_id] = []
                for o in o_list:
                    msg_dict[record_type][main_id].append(o)
                    if "replacement_id_list" in o:
                        replacor_id = o["replacement_id_list"][0]
                        msg = "inherits %s as of data release %s" % (main_id, rel)
                        oo = {"type":"inherits", "inherited_id":main_id,"description":msg}
                        if replacor_id not in msg_dict[record_type]:
                            msg_dict[record_type][replacor_id] = []
                        msg_dict[record_type][replacor_id].append(oo)

    record_count = 0
    record_count += output_glygen(msg_dict)
    record_count += output_glytoucan_only(history_dict, glygen_ac_dict, glytoucan_ac_dict)
    record_count += output_uniprotkb_only(history_dict, glygen_ac_dict,ac2canon)
    record_count += output_batches(batch_doc_list)

    out_json = {}
    for rel in release_list:
        for record_type in ["protein", "glycan"]:
            for tax_id in sorted(ref_tax_id_list) + ["unknown"]:
                tax_name = species_obj[tax_id]["common_name"] if tax_id in species_obj else tax_id
                added, discontinued_in_glygen = 0, 0
                if rel in record_stat:
                    if record_type in record_stat[rel]:
                        if tax_id in record_stat[rel][record_type]:
                            added = record_stat[rel][record_type][tax_id]["added"]
                            discontinued_in_glygen = record_stat[rel][record_type][tax_id]["discontinued_in_glygen"]
                o = {"release":rel, "release_date":ts_dict[rel], "added":added, "discontinued_in_glygen":discontinued_in_glygen}
                if tax_name not in out_json:
                    out_json[tax_name] = {}
                if record_type not in out_json[tax_name]:
                    out_json[tax_name][record_type] = []
                out_json[tax_name][record_type].append(o)


    msg = "make-idtrackdb: final created: %s idtrack objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "a")




if __name__ == '__main__':
    main()

