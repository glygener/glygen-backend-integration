#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob


def get_protein_sec_stats(doc, record_type):

    tmp_dict = {
        "protein":{
            "glycosylation":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites"
            ],
            "glycosylation_reported":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites"
            ]
            ,"glycosylation_reported_with_glycans":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites"
            ]
            ,"glycosylation_predicted":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites"
            ]
            ,"glycosylation_automatic_literature_mining":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites"
            ]
            ,"phosphorylation":["total", "total_sites"]
            ,"snv":["total", "total_sites"]
            ,"snv_disease":["total", "total_sites"]
            ,"snv_non_disease":["total", "total_sites"]
            ,"publication":["total"]
        },
        "site":{
            "glycosylation":["total", "n_linked_glycans",  "o_linked_glycans"]
            ,"glycosylation_reported":["total", "n_linked_glycans", "o_linked_glycans"]
            ,"glycosylation_reported_with_glycans":["total", "n_linked_glycans", "o_linked_glycans"]
            ,"glycosylation_predicted":["total", "n_linked_glycans", "o_linked_glycans"]
            ,"glycosylation_automatic_literature_mining":["total", "n_linked_glycans", "o_linked_glycans"]
            ,"phosphorylation":["total"]
            ,"snv":["total"]
            ,"snv_disease":["total"]
            ,"snv_non_disease":["total"]
            ,"publication":["total"]
        }
    }
    
    stats_obj = {}
    for table_id in tmp_dict[record_type]:
        obj = {}
        for f in tmp_dict[record_type][table_id]:
            obj[f] = 0    
        stats_obj[table_id] = obj


    for table_id in ["phosphorylation", "publication"]:
        if table_id not in doc:
            continue
        if doc[table_id] != []:
            for obj in doc[table_id]:
                stats_obj[table_id]["total"] += 1
                if "start_pos" in obj:
                    stats_obj[table_id]["total_sites"] += 1
   
    table_id_all = "snv" 
    if table_id_all in doc:
        for obj in doc[table_id_all]:
            stats_obj[table_id_all]["total"] += 1
            if "start_pos" in obj:
                stats_obj[table_id_all]["total_sites"] += 1
            if "disease" in obj["keywords"]:
                table_id = "snv_disease"
                stats_obj[table_id]["total"] += 1
                if "start_pos" in obj:
                    stats_obj[table_id]["total_sites"] += 1 
            else:
                table_id = "snv_non_disease"
                stats_obj[table_id]["total"] += 1
                if "start_pos" in obj:
                    stats_obj[table_id]["total_sites"] += 1 
 
    if "glycosylation" in doc:
        if doc["glycosylation"] != []:
            seen = {}
            for obj in doc["glycosylation"]:
                site_lbl = obj["site_lbl"] if "site_lbl" in obj else ""
                glytoucan_ac = obj["saccharide"] if "saccharide" in obj else ""
                if obj["site_category"] == "reported":
                    table_id = "glycosylation_reported"
                    table_id_all = table_id.split("_")[0]
                    stats_obj[table_id]["total"] += 1
                    stats_obj[table_id_all]["total"] += 1
                    if site_lbl not in seen:
                        stats_obj[table_id]["total_sites"] += 1
                        stats_obj[table_id_all]["total_sites"] += 1
 
                    #is_site = False
                    #if "start_pos" in obj and record_type == "protein":
                    #    stats_obj[table_id]["total_sites"] += 1
                    #    stats_obj[table_id_all]["total_sites"] += 1
                    #    is_site = True
                    if glytoucan_ac != "" and glytoucan_ac not in seen and obj["type"].lower() == "n-linked":
                        stats_obj[table_id]["n_linked_glycans"] += 1
                        stats_obj[table_id_all]["n_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["n_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["n_linked_glycan_sites"] += 1
                    if obj["type"].lower() == "o-linked":
                        stats_obj[table_id]["o_linked_glycans"] += 1
                        stats_obj[table_id_all]["o_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["o_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["o_linked_glycan_sites"] += 1
                elif obj["site_category"] == "reported_with_glycan":
                    table_id = "glycosylation_reported_with_glycans"
                    table_id_all = table_id.split("_")[0]
                    stats_obj[table_id]["total"] += 1
                    stats_obj[table_id_all]["total"] += 1
                    is_site = False
                    if "start_pos" in obj and record_type == "protein":
                        stats_obj[table_id]["total_sites"] += 1
                        stats_obj[table_id_all]["total_sites"] += 1
                        is_site = True
                    if obj["type"].lower() == "n-linked":
                        stats_obj[table_id]["n_linked_glycans"] += 1
                        stats_obj[table_id_all]["n_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["n_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["n_linked_glycan_sites"] += 1
                    if obj["type"].lower() == "o-linked":
                        stats_obj[table_id]["o_linked_glycans"] += 1
                        stats_obj[table_id_all]["o_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["o_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["o_linked_glycan_sites"] += 1
                elif obj["site_category"] == "automatic_literature_mining":
                    table_id = "glycosylation_automatic_literature_mining"
                    table_id_all = table_id.split("_")[0]
                    stats_obj[table_id]["total"] += 1
                    stats_obj[table_id_all]["total"] += 1
                    is_site = False
                    if "start_pos" in obj and record_type == "protein":
                        stats_obj[table_id]["total_sites"] += 1
                        stats_obj[table_id_all]["total_sites"] += 1
                        is_site = True
                    if obj["type"].lower() == "n-linked":
                        stats_obj[table_id]["n_linked_glycans"] += 1
                        stats_obj[table_id_all]["n_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["n_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["n_linked_glycan_sites"] += 1
                    if obj["type"].lower() == "o-linked":
                        stats_obj[table_id]["o_linked_glycans"] += 1
                        stats_obj[table_id_all]["o_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["o_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["o_linked_glycan_sites"] += 1
                elif obj["site_category"] == "predicted":
                    table_id = "glycosylation_predicted"
                    table_id_all = table_id.split("_")[0]
                    stats_obj[table_id]["total"] += 1
                    stats_obj[table_id_all]["total"] += 1
                    is_site = False
                    if "start_pos" in obj and record_type == "protein":
                        stats_obj[table_id]["total_sites"] += 1
                        stats_obj[table_id_all]["total_sites"] += 1
                        is_site = True
                    if obj["type"].lower() == "n-linked":
                        stats_obj[table_id]["n_linked_glycans"] += 1
                        stats_obj[table_id_all]["n_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["n_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["n_linked_glycan_sites"] += 1
                    if obj["type"].lower() == "o-linked":
                        stats_obj[table_id]["o_linked_glycans"] += 1
                        stats_obj[table_id_all]["o_linked_glycans"] += 1
                        if is_site == True and record_type == "protein":
                            stats_obj[table_id]["o_linked_glycan_sites"] += 1
                            stats_obj[table_id_all]["o_linked_glycan_sites"] += 1
                seen[site_lbl] = True

    #doc["section_stats"] = stats_obj
    doc["section_stats"] = []
    for table_id in stats_obj:
        obj = {"table_id":table_id, "table_stats":[]}
        for f in stats_obj[table_id]:
            obj["table_stats"].append({"field":f, "count":stats_obj[table_id][f]})
        doc["section_stats"].append(obj)
    
 
    return








def get_sec_stats(doc, record_type):

    sec_map = {
        "glycan":{
            "glycoprotein":["total", "total_sites"]
            ,"expression":["total", "total_sites"]
            ,"expression_tissue":["total", "total_sites"]
            ,"expression_cell_line":["total", "total_sites"]
            ,"publication":["total"]
        }
        ,"motif":{
            "glycans":["total"]
        }
        ,"site":{
            "glycosylation":["total"]
            ,"phosphorylation":["total"]
            ,"snv":["total"]
        }
        ,"biomarker":{
            "disease":["total"]
            ,"biomarker_components":["total"]
            ,"publication":["total"]
        }
        ,"publication":{
            "referenced_proteins":["total"]
            ,"referenced_glycans":["total"]
            ,"glycosylation":[]
            ,"phosphorylation":[]
        }
    }
    tmp_dict = sec_map[record_type]
    stats_obj = {}
    for table_id in tmp_dict:
        obj = {}
        for f in tmp_dict[table_id]:
            obj[f] = 0    
        stats_obj[table_id] = obj

    for table_id in tmp_dict:
        if record_type == "publication" and table_id in ["glycosylation", "phosphorylation"]:
            tmp_doc = {table_id:doc[table_id]}
            get_protein_sec_stats(tmp_doc, "protein")
            stats_obj[table_id] = {}
            for obj in tmp_doc["section_stats"]:
                if obj["table_id"].split("_")[0] != table_id:
                    continue
                stats_obj[obj["table_id"]] = {}
                for o in obj["table_stats"]:
                    stats_obj[obj["table_id"]][o["field"]] = o["count"]
        elif record_type == "glycan" and table_id == "expression_tissue":
            table_id_all = table_id.split("_")[0]
            for obj in doc[table_id_all]:
                if obj["tissue"]  != {}:
                    stats_obj[table_id]["total"] += 1
                    if "start_pos" in obj:
                        stats_obj[table_id]["total_sites"] += 1
        elif record_type == "glycan" and table_id == "expression_cell_line":
            table_id_all = table_id.split("_")[0]
            for obj in doc[table_id_all]:
                if obj["cell_line"]  != {}: 
                    stats_obj[table_id]["total"] += 1
                    if "start_pos" in obj:
                        stats_obj[table_id]["total_sites"] += 1
        elif doc[table_id] != []:
            for obj in doc[table_id]:
                stats_obj[table_id]["total"] += 1
                if "start_pos" in obj:
                    stats_obj[table_id]["total_sites"] += 1
    
    #doc["section_stats"] = stats_obj
    doc["section_stats"] = []
    for table_id in stats_obj:
        obj = {"table_id":table_id, "table_stats":[]}
        if record_type == "motif" and table_id == "glycans":
            obj["table_id"] = "results"
        for f in stats_obj[table_id]:
            obj["table_stats"].append({"field":f, "count":stats_obj[table_id][f]})
        doc["section_stats"].append(obj)
    
    
    return


