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
import csvutil



def get_children_tree(p_id):

    obj_list = [] 
    if p_id in tree_dict:
        for c_id in tree_dict[p_id]:
            name = name_dict[c_id] if c_id in name_dict else ""
            o = {"id":c_id, "label":name, "children":[]}
            o["children"] = get_children_tree(c_id)
            #if o["children"] == []:
            #    o.pop("children")
            obj_list.append(o)
         

    return obj_list



def get_flat_dict(p_id):
        
    name = name_dict[p_id] if p_id in name_dict else ""
    o = {"id":p_id, "label":name}
    obj_list = [o]
    if p_id in tree_dict:
        for c_id in tree_dict[p_id]:
            name = name_dict[c_id] if c_id in name_dict else ""
            o = {"id":c_id, "label":name}
            obj_list.append(o)
            obj_list += get_flat_dict(c_id)

    return obj_list


def load_tree_dict():

    name_dict = {}
    tree_dict = {}
    data_frame = {}
    in_file = "reviewed/protein_disease_tree.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        c_id = row[f_list.index("disease_id")].replace("DOID_", "DOID:")
        c_name = row[f_list.index("disease_name")]
        p_id = row[f_list.index("parent_disease_id")].replace("DOID_", "DOID:")
        p_name = row[f_list.index("parent_disease_name")]
        name_dict[c_id] = c_name
        name_dict[p_id] = p_name
        if p_id not in tree_dict:
            tree_dict[p_id] = {}
        tree_dict[p_id][c_id] = True


    return tree_dict, name_dict




def load_disease_idmap():

    map_dict = {}
    csvutil.load_dictionaries(map_dict, "generated/misc/")
    


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
            xref_id = do_id.split(":")[1]
            xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
            xref_badge = libgly.get_xref_badge(map_dict,xref_key)
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
                        xref_id = mapped_id.split(":")[1]
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                        xref_badge = libgly.get_xref_badge(map_dict,xref_key)
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
                        xref_id = mapped_id.split(":")[1]
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                        xref_badge = libgly.get_xref_badge(map_dict,xref_key)
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
                xref_id = mondo_id.split(":")[1]
                xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict,xref_key)
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
                        xref_id = mapped_id.split(":")[1]
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                        xref_badge = libgly.get_xref_badge(map_dict,xref_key)
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
                xref_id = mim_id.split(":")[1]
                xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict,xref_key)
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
    global main_dict
    global tree_dict
    global name_dict


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    path_obj["reviewed"] = "reviewed/"

    
    global is_cited



    tree_dict, name_dict = load_tree_dict()

    flat_dict = {}    
    children_dict = {}
    for p_id in tree_dict:
        children_dict[p_id] = get_children_tree(p_id)
        flat_dict[p_id] = get_flat_dict(p_id)



    is_cited = libgly.get_is_cited()


    data_dir = "reviewed/"

    final_dict = load_disease_idmap()
    record_count = 0
    main_id_list = list(final_dict.keys())
    for main_id in main_id_list:
        name = final_dict[main_id]["recommended_name"]["name"]
        final_dict[main_id]["children"] = children_dict[main_id] if main_id in children_dict else []
        id_list, name_list = [main_id], [name]
        if main_id in flat_dict:
            for o in flat_dict[main_id]:
                id_list.append(o["id"])
                name_list.append(o["label"])
        
        final_dict[main_id]["id_list"] = list(set(id_list))
        final_dict[main_id]["name_list"]  = list(set(name_list))
        out_file = "jsondb/interm/diseasedb/%s.json" % (main_id.replace(":",".").lower())
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(final_dict[main_id], indent=4)))
        record_count += 1


    log_file = "logs/make-diseasedb.log"
    msg = "make-diseasedb: final created: %s disease objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "w")



if __name__ == '__main__':
    main()

