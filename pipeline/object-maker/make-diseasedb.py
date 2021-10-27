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


def load_dictionaries(map_dict, misc_dir):

    dict_list_obj = json.loads(open("../../conf/protein_dictionaries.json", "r").read())
    for dict_name in dict_list_obj:
        map_dict[dict_name] = {}
        ind_list = dict_list_obj[dict_name]["indexlist"]
        for pattern in dict_list_obj[dict_name]["fileglob"]:
            for in_file in glob.glob(misc_dir + pattern):
                sheet_obj = {}
                libgly.load_sheet(sheet_obj, in_file, ",")
                for row in sheet_obj["data"]:
                    if row ==[] or row[ind_list[0]][0] == "#":
                        continue
                    key = row[ind_list[0]]
                    val = row[ind_list[1]]
                    if key not in map_dict[dict_name]:
                        map_dict[dict_name][key] = []
                    map_dict[dict_name][key].append(val)

    return


def load_disease_idmap():

    map_dict = {}
    load_dictionaries(map_dict, "generated/misc/")
    



    master_dict = {}
    file_list = glob.glob(path_obj["reviewed"] + "/*_protein_disease_*.csv")
    file_list += glob.glob(path_obj["reviewed"] + "/*_protein_expression_disease.csv")
    file_list += glob.glob(path_obj["reviewed"] + "/*_protein_mutation_cance.csv")
    for in_file in file_list:
        if in_file.find(".stat.csv") != -1:
            continue
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        lcount = 1
        for row in data_frame["data"]:
            do_id, mondo_id, mim_id = "", "", ""
            do_id = row[f_list.index("do_id")] if "do_id" in f_list else do_id
            mondo_id = row[f_list.index("mondo_id")] if "mondo_id" in f_list else mondo_id
            mim_id = row[f_list.index("mim_id")] if "mim_id" in f_list else mim_id
            do_id = do_id.replace("DOID:", "")
            mondo_id = mondo_id.replace("MONDO:", "")
            if do_id != "":
                main_id = "DOID:%s" % (do_id)
                master_dict[main_id] = True
            if mondo_id != "":
                main_id = "MONDO:%s" % (mondo_id)
                master_dict[main_id] = True
            if mim_id != "":
                main_id = "MIM:%s" % (mim_id)
                master_dict[main_id] = True




    disease_idmap = {}
    is_mapped = {}
    data_frame = {}
    in_file = path_obj["reviewed"]+ "/protein_disease_idmap.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        mondo_id = row[f_list.index("mondo_id")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        mim_id = xref_id if xref_key == "protein_xref_omim" else ""
        if do_id != "":
            main_id = "DOID:%s" % (do_id)
            if main_id not in disease_idmap:
                disease_idmap[main_id] = {"do2mondo":[], "do2mim":[]}
            if mondo_id != "" and mondo_id not in disease_idmap[main_id]["do2mondo"]:
                m_id = "MONDO:%s"%(mondo_id)
                disease_idmap[main_id]["do2mondo"].append(m_id)
                is_mapped[m_id] = True 
            if mim_id != "" and mim_id not in disease_idmap[main_id]["do2mim"]:
                m_id = "MIM:%s"%(mim_id)
                disease_idmap[main_id]["do2mim"].append(m_id)
                is_mapped[m_id] = True
        elif mondo_id != "":
            main_id = "MONDO:%s" % (mondo_id)
            if main_id not in disease_idmap:
                disease_idmap[main_id] = {"mondo2mim":[]}
            if mim_id != "" and mim_id not in disease_idmap[main_id]["mondo2mim"]:
                m_id = "MIM:%s"%(mim_id)
                disease_idmap[main_id]["mondo2mim"].append(m_id)
                is_mapped[m_id] = True


    #add unmapped IDs
    for main_id in master_dict:
        if main_id not in disease_idmap:
            disease_idmap[main_id] = {"mondo2mim":[], "do2mim":[], "do2mondo":[]}



    rec_name_dict = {}
    syn_name_dict = {}
    data_frame = {}
    in_file = path_obj["reviewed"] + "/protein_disease_names.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    tmp_fl = data_frame["fields"]
    for tmp_row in data_frame["data"]:
        xref_key = tmp_row[tmp_fl.index("xref_key")]
        xref_id = tmp_row[tmp_fl.index("xref_id")]
        name_type = tmp_row[tmp_fl.index("name_type")]
        name = tmp_row[tmp_fl.index("name")]
        disease_desc = tmp_row[tmp_fl.index("description")]
        pref = xref_key.upper()
        pref = "DOID" if pref == "DO" else pref
        pref = "MIM" if pref == "OMIM" else pref
        main_id = "%s:%s" % (pref, xref_id)
        if name_type == "recommended_name":
            if main_id not in rec_name_dict:
                rec_name_dict[main_id] = {}
            rec_name_dict[main_id][name] = disease_desc
        else:
            if main_id not in syn_name_dict:
                syn_name_dict[main_id] = {}
            syn_name_dict[main_id][name] = disease_desc


    final_dict = {}
    for main_id in disease_idmap:
        if main_id in is_mapped:
            continue
        if main_id.find("DOID:") != -1:
            do_id = main_id
            xref_key = "protein_xref_do"
            xref_url = map_dict["xrefkey2url"][xref_key][0] % (do_id.split(":")[1])
            xref_badge = map_dict["xrefkey2badge"][xref_key][0]
            if do_id in rec_name_dict:
                final_dict[main_id] = {"disease_id":do_id, "recommended_name":{}, "synonyms": []}
                for name in rec_name_dict[do_id]:
                    desc = rec_name_dict[do_id][name]
                    o = {"id":do_id, "resource":xref_badge,"url":xref_url,
                            "name":name, "description":desc}
                    final_dict[main_id]["recommended_name"] = o
                #add DO synonyms
                if do_id in syn_name_dict:
                    for name in syn_name_dict[do_id]:
                        desc = syn_name_dict[do_id][name]
                        o = {"id":do_id, "resource":xref_badge,"url":xref_url, 
                                "name":name, "description":desc}
                        final_dict[main_id]["synonyms"].append(o)
                
                #add MONDO synonyms
                for mapped_id in disease_idmap[do_id]["do2mondo"]:
                    if mapped_id in rec_name_dict:
                        xref_key = "protein_xref_mondo"
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (mapped_id.split(":")[1])
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                        for name in rec_name_dict[mapped_id]:
                            desc = rec_name_dict[mapped_id][name]
                            o = {"id":mapped_id, "resource":xref_badge,"url":xref_url,  
                                    "name":name, "description":desc}
                            final_dict[main_id]["synonyms"].append(o)
                    if mapped_id in syn_name_dict:
                        for name in syn_name_dict[mapped_id]:
                            desc = syn_name_dict[mapped_id][name]
                            o = {"id":mapped_id, "resource":xref_badge,"url":xref_url,
                                    "name":name, "description":desc}
                            final_dict[main_id]["synonyms"].append(o)
                #add MIM synonyms
                for mapped_id in disease_idmap[do_id]["do2mim"]:
                    if mapped_id in rec_name_dict:
                        xref_key = "protein_xref_omim"
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (mapped_id.split(":")[1])
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                        for name in rec_name_dict[mapped_id]:
                            desc = rec_name_dict[mapped_id][name]
                            o = {"id":mapped_id, "resource":xref_badge,"url":xref_url,
                                    "name":name, "description":desc}
                            final_dict[main_id]["synonyms"].append(o)
                    if mapped_id in syn_name_dict:
                        for name in syn_name_dict[mapped_id]:
                            desc = syn_name_dict[mapped_id][name]
                            o = {"id":mapped_id, "resource":xref_badge,"url":xref_url,
                                    "name":name, "description":desc}
                            final_dict[main_id]["synonyms"].append(o)


        if main_id.find("MONDO:") != -1:
            mondo_id = main_id
            if mondo_id in rec_name_dict:
                xref_key = "protein_xref_mondo"
                xref_url = map_dict["xrefkey2url"][xref_key][0] % (mondo_id.split(":")[1])
                xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                final_dict[main_id] = {"disease_id":mondo_id, "recommended_name":{}, "synonyms": []}
                for name in rec_name_dict[mondo_id]:
                    desc = rec_name_dict[mondo_id][name]
                    o = {"id":mondo_id, "resource":xref_badge,"url":xref_url,
                            "name":name, "description":desc}
                    final_dict[main_id]["recommended_name"] = o
                if mondo_id in syn_name_dict:
                    for name in syn_name_dict[mondo_id]:
                        desc = syn_name_dict[mondo_id][name]
                        o = {"id":mondo_id, "resource":xref_badge,"url":xref_url,
                                "name":name, "description":desc}
                        final_dict[main_id]["synonyms"].append(o)

                for mapped_id in disease_idmap[mondo_id]["mondo2mim"]:
                    if mapped_id in rec_name_dict:
                        xref_key = "protein_xref_omim"
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (mapped_id.split(":")[1])
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                        for name in rec_name_dict[mapped_id]:
                            desc = rec_name_dict[mapped_id][name]
                            o = {"id":mapped_id, "resource":xref_badge,"url":xref_url,
                                  "name":name, "description":desc}
                            final_dict[main_id]["synonyms"].append(o)
                    if mapped_id in syn_name_dict:
                        for name in syn_name_dict[mapped_id]:
                            desc = syn_name_dict[mapped_id][name]
                            o = {"id":mapped_id, "resource":xref_badge,"url":xref_url,
                                 "name":name, "description":desc}
                            final_dict[main_id]["synonyms"].append(o)
         
        if main_id.find("MIM:") != -1:
            mim_id = main_id
            if mim_id in rec_name_dict:
                xref_key = "protein_xref_omim"
                xref_url = map_dict["xrefkey2url"][xref_key][0] % (mim_id.split(":")[1])
                xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                final_dict[main_id] = {"disease_id":mim_id, "recommended_name":{}, "synonyms": []}
                for name in rec_name_dict[mim_id]:
                    desc = rec_name_dict[mim_id][name]
                    o = {"id":mim_id, "resource":xref_badge,"url":xref_url,
                            "name":name, "description":desc}
                    final_dict[main_id]["recommended_name"] = o
                if mim_id in syn_name_dict:
                    for name in syn_name_dict[mim_id]:
                        desc = syn_name_dict[mim_id][name]
                        o = {"id":mim_id, "resource":xref_badge,"url":xref_url,
                                "name":name, "description":desc}
                        final_dict[main_id]["synonyms"].append(o)



    return final_dict



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

    final_dict = load_disease_idmap()
    record_count = 0
    for main_id in final_dict:
        out_file = "jsondb/diseasedb/%s.json" % (main_id.replace(":",".").lower())
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(final_dict[main_id], indent=4)))
        record_count += 1

    print "make-diseasedb: final created: %s disease objects" % (record_count)




if __name__ == '__main__':
    main()

