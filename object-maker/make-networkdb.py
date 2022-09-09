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


def add_decendant(tree_dict, src_type, src_id, dst_type, dst_id):

    if src_type not in tree_dict:
        tree_dict[src_type] = {}
    if src_id not in tree_dict[src_type]:
        tree_dict[src_type][src_id] = {}
    if dst_type not in tree_dict[src_type][src_id]:
        tree_dict[src_type][src_id][dst_type] = {}
    tree_dict[src_type][src_id][dst_type][dst_id] = True

    if dst_type not in tree_dict:
        tree_dict[dst_type] = {}
    if dst_id not in tree_dict[dst_type]:
        tree_dict[dst_type][dst_id] = {}
    if src_type not in tree_dict[dst_type][dst_id]:
        tree_dict[dst_type][dst_id][src_type] = {}
    tree_dict[dst_type][dst_id][src_type][src_id] = True


    return





def get_decendants(subtree_obj, src_id):

    out_obj = {}
    for dec_type in subtree_obj[src_id]:
        if dec_type not in out_obj:
            out_obj[dec_type] = {}
        for dec_id in subtree_obj[src_id][dec_type]:
            out_obj[dec_type][dec_id] = True
        for child_id in subtree_obj[src_id][dec_type]:
            if child_id in subtree_obj:
                tmp_obj = get_decendants(subtree_obj, child_id)
                for child_dec_type in tmp_obj:
                    for child_dec_id in tmp_obj[child_dec_type]:
                        if child_dec_type not in out_obj:
                            out_obj[child_dec_type] = {}
                        out_obj[child_dec_type][child_dec_id] = True


    return out_obj


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

    
    #glycan_file_list = glob.glob("jsondb/glycandb/G1768*.json")
    #protein_file_list = glob.glob("jsondb/proteindb/P1421*.json")
    glycan_file_list = glob.glob("jsondb/glycandb/*.json")
    protein_file_list = glob.glob("jsondb/proteindb/*.json")
            


    seen_edge = {}
    tree_dict = {}

    record_count = 0
    for json_file in glycan_file_list:
        doc = json.loads(open(json_file,"r").read())
        src_type, src_id = "glycan", doc["glytoucan_ac"]
        linked_flag = False
        sec = "species"
        if doc[sec] != []:
            linked_flag = True
            for obj in doc[sec]:
                dst_type, dst_id = "species", obj["taxid"]
                add_decendant(tree_dict, src_type, src_id, dst_type, dst_id)
        sec = "motifs"
        if doc[sec] != []:
            linked_flag = True
            for obj in doc[sec]:
                dst_type, dst_id = "motif", obj["id"]
                add_decendant(tree_dict, src_type, src_id, dst_type, dst_id)
        
        sec = "enzyme"
        if doc[sec] != []:
            linked_flag = True
            for obj in doc[sec]:
                canon = obj["uniprot_canonical_ac"]
                dst_type, dst_id = "enzyme", "enzyme." + canon
                add_decendant(tree_dict, src_type, src_id, dst_type, dst_id)


        if linked_flag == False:
            if src_type not in tree_dict:
                tree_dict[src_type] = {}
            tree_dict[src_type][src_id] = {}

        record_count += 1



    record_count = 0
    for json_file in protein_file_list:
        doc = json.loads(open(json_file,"r").read())
        canon = doc["uniprot_canonical_ac"]
        tax_id = doc["species"][0]["taxid"]
        #gene_name = doc["gene"][0]["name"]
        gene_name = "gene." + canon

        
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
                    site_id = "%s.%s.%s" % (canon,start_pos, end_pos)
                    src_type, src_id = "protein",canon
                    dst_type, dst_id = "site", site_id
                    add_decendant(tree_dict, src_type, src_id, dst_type, dst_id)
                    
                    #src_type, src_id = "site",site_id
                    #dst_type, dst_id = "species", tax_id
                    #add_decendant(tree_dict, src_type, src_id, dst_type, dst_id)
                    
                    if sec == "glycosylation":
                        if obj["site_category"] =="reported_with_glycan":
                            src_type, src_id = "site", site_id
                            dst_type, dst_id = "glycan", obj["glytoucan_ac"]
                            add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)

        sec = "interactions"
        if doc[sec] != []:
            for obj in doc[sec]:
                if obj["interaction_type"].lower() == "glycosaminoglycan":
                    src_type, src_id = "protein",canon
                    dst_type, dst_id = "glycan", obj["interactor_id"]
                    add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)
        

        sec = "synthesized_glycans"
        if doc[sec] != []:
            for obj in doc[sec]:
                src_type, src_id = "protein", canon
                dst_type, dst_id = "glycan", obj["glytoucan_ac"]
                add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)

        src_type, src_id = "protein", canon
        dst_type, dst_id = "species", tax_id
        add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)
    

        src_type, src_id = "protein", canon
        dst_type, dst_id = "gene", gene_name
        add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)

        src_type, src_id = "gene", gene_name
        dst_type, dst_id = "species", tax_id
        add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)



        sec = "keywords"
        if doc[sec] != []:
            if "enzyme" in doc[sec]:
                src_type, src_id = "protein", canon
                dst_type, dst_id = "enzyme", "enzyme." + canon
                add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)

        #diseases from snv are linked both to protein and site        
        sec = "snv"
        for obj in doc[sec]:
            if "disease" not in obj:
                continue
            site_id = "%s.%s.%s" % (canon,obj["start_pos"], obj["end_pos"])

            for o in obj["disease"]:
                src_type, src_id = "protein", canon
                dst_type, dst_id = "disease",o["recommended_name"]["id"] 
                add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)
                src_type, src_id = "gene", gene_name
                dst_type, dst_id = "disease",o["recommended_name"]["id"]
                add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)

                src_type, src_id = "site",site_id
                dst_type, dst_id = "disease",o["recommended_name"]["id"] 
                add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)
                

        #diseases from expression_disease are linked to protein
        sec = "expression_disease"
        for obj in doc[sec]:
            if "disease" not in obj:
                continue
            for o in obj["disease"]:
                src_type, src_id = "protein", canon
                dst_type, dst_id = "disease",o["recommended_name"]["id"] 
                add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)
                src_type, src_id = "gene", gene_name
                dst_type, dst_id = "disease",o["recommended_name"]["id"]
                add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)


        #diseases from disease are linked to protein
        sec = "disease"
        for obj in doc[sec]:
            src_type, src_id = "protein", canon
            dst_type, dst_id = "disease",obj["recommended_name"]["id"]
            add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)
       
            src_type, src_id = "gene", gene_name
            dst_type, dst_id = "disease",obj["recommended_name"]["id"]
            add_decendant(tree_dict, src_type,src_id, dst_type, dst_id)

        record_count += 1



    
    batch_count = 0
    record_count = 0
    obj_list = []
    for src_type in tree_dict:
        batch_size = 1 if src_type == "species" else 50
        for src_id in tree_dict[src_type]:
            if record_count > 0 and record_count%batch_size == 0:
                batch_count += 1
                doc = {"batch":batch_count,"recordlist":obj_list}
                out_file = "jsondb/networkdb/batch.%s.json" % (batch_count)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(doc, indent=4)))
                obj_list = []
            obj = {"record_id":src_id, "record_type":src_type, "linkeage":{}}
            seen_dec = get_decendants(tree_dict[src_type], src_id)
            for dec_type in seen_dec:
                obj["linkeage"][dec_type] = list(seen_dec[dec_type].keys())
            obj_list.append(obj)
            record_count += 1

    batch_count += 1
    doc = {"batch":batch_count,"recordlist":obj_list}
    out_file = "jsondb/networkdb/batch.%s.json" % (batch_count)
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(doc, indent=4)))

    print ("make-networkdb: ... final created: %s network objects" % (record_count))


if __name__ == '__main__':
    main()
