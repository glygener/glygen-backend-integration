#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob


def get_sort_fields(obj, tmp_list):
    if type(obj) is dict:
        for f in obj:
            if type(obj[f]) in [int, str, float]:
                if f not in tmp_list:
                    tmp_list.append(f)
            elif type(obj[f]) is dict:
                for ff in obj[f]:
                    if type(obj[f][ff]) in [int, str, float]:
                        cmb = f + "." + ff
                        if cmb not in tmp_list:
                            tmp_list.append(cmb)

    return


def get_protein_sec_stats(doc, record_type):

    tmp_dict = {
        "protein":{
            "glycosylation":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites", "n_linked_sites", "o_linked_sites" 
            ],
            "glycosylation_reported":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites", "n_linked_sites", "o_linked_sites"
            ]
            ,"glycosylation_reported_with_glycan":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites", "n_linked_sites", "o_linked_sites"
            ]
            ,"glycosylation_predicted":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites", "n_linked_sites", "o_linked_sites"
            ]
            ,"glycosylation_automatic_literature_mining":[
                "total", "total_sites", "n_linked_glycans", "n_linked_glycan_sites", "o_linked_glycans",
                "o_linked_glycan_sites", "n_linked_sites", "o_linked_sites"
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
            ,"glycosylation_reported_with_glycan":["total", "n_linked_glycans", "o_linked_glycans"]
            ,"glycosylation_predicted":["total", "n_linked_glycans", "o_linked_glycans"]
            ,"glycosylation_automatic_literature_mining":["total", "n_linked_glycans", "o_linked_glycans"]
            ,"phosphorylation":["total"]
            ,"snv":["total"]
            ,"snv_disease":["total"]
            ,"snv_non_disease":["total"]
            ,"publication":["total"]
        }
    }
   
    sort_obj = {} 
    stats_obj = {}
    for table_id in tmp_dict[record_type]:
        obj = {}
        for f in tmp_dict[record_type][table_id]:
            obj[f] = 0    
        stats_obj[table_id] = obj
        sort_obj[table_id] = []
    #print (json.dumps(stats_obj, indent=4))
    #exit()


    for table_id_all in ["phosphorylation", "publication"]:
        if table_id_all not in doc:
            continue
        if doc[table_id_all] != []:
            seen = {}
            for obj in doc[table_id_all]:
                #sort_obj[table_id_all] = get_sort_fields(obj, sort_obj[table_id_all])
                get_sort_fields(obj, sort_obj[table_id_all])
                site_lbl_all = table_id_all + "_" + obj["site_lbl"] if "site_lbl" in obj else ""
                stats_obj[table_id_all]["total"] += 1
                if "total_sites" in stats_obj[table_id_all]:
                    if site_lbl_all not in seen:
                        stats_obj[table_id_all]["total_sites"] += 1
                        seen[site_lbl_all] = True

    table_id_all = "snv" 
    if table_id_all in doc:
        seen = {}
        for obj in doc[table_id_all]:
            #sort_obj[table_id_all] = get_sort_fields(obj, sort_obj[table_id_all])
            get_sort_fields(obj, sort_obj[table_id_all])
            site_lbl_all = table_id_all + "_" + obj["site_lbl"] if "site_lbl" in obj else ""
            stats_obj[table_id_all]["total"] += 1
            if site_lbl_all not in seen and "total_sites" in stats_obj[table_id_all]:
                stats_obj[table_id_all]["total_sites"] += 1
                seen[site_lbl_all] = True
            if "disease" in obj["keywords"]:
                table_id = "snv_disease"
                sort_obj[table_id] = sort_obj[table_id_all] 
                site_lbl = table_id + "_" + obj["site_lbl"] if "site_lbl" in obj else ""
                stats_obj[table_id]["total"] += 1
                if site_lbl not in seen and "total_sites" in stats_obj[table_id]:
                    stats_obj[table_id]["total_sites"] += 1 
                    seen[site_lbl] = True
            else:
                table_id = "snv_non_disease"
                sort_obj[table_id] = sort_obj[table_id_all] 
                site_lbl = table_id + "_" + obj["site_lbl"] if "site_lbl" in obj else ""
                stats_obj[table_id]["total"] += 1
                if site_lbl not in seen and "total_sites" in stats_obj[table_id]:
                    stats_obj[table_id]["total_sites"] += 1 
                    seen[site_lbl] = True

    total_dict = {}
    cmb_dict = {} 
    seen_gtype = {}
    if "glycosylation" in doc:
        if doc["glycosylation"] != []:
            for obj in doc["glycosylation"]:
                table_id_all = "glycosylation"
                if table_id_all not in total_dict:
                    total_dict[table_id_all] = 0
                total_dict[table_id_all] += 1
                get_sort_fields(obj, sort_obj[table_id_all])
                glytoucan_ac = obj["glytoucan_ac"] if "glytoucan_ac" in obj else ""
                is_n_linked = obj["type"].lower() == "n-linked"
                is_o_linked = obj["type"].lower() == "o-linked"
                gtype = obj["type"].lower()
                seen_gtype[gtype] = True
                start_pos = obj["start_pos"] if "start_pos" in obj else "*"
                site_type = "fuzzy" if start_pos == "*" else "exact"
                update_cmb_dict(cmb_dict, table_id_all, site_type, start_pos, glytoucan_ac, gtype)
                for site_cat in obj["site_category_dict"]:
                    if site_cat in ["reported_with_glycan", "reported", "automatic_literature_mining", "predicted"]:
                        table_id = "glycosylation_" + site_cat
                        update_cmb_dict(cmb_dict, table_id, site_type, start_pos, glytoucan_ac, gtype)
                        sort_obj[table_id] = sort_obj[table_id_all] 
                        if table_id not in total_dict:
                            total_dict[table_id] = 0
                        total_dict[table_id] += 1 
    #print (json.dumps(cmb_dict, indent=4))
    
    for table_id in cmb_dict:
        if cmb_dict[table_id]["exact"]["one"] != {}:
            stats_obj[table_id] = {}
            #stats_obj[table_id]["exact_sites"] = len(cmb_dict[table_id]["exact"]["one"].keys())
            for s in ["exact","fuzzy"]:
                #s_lbl = s if s == "fuzzy_" else "xxx"
                s_lbl = s + "_"
                k = "%ssites"%(s_lbl)
                stats_obj[table_id][k] = len(cmb_dict[table_id][s]["one"].keys())
                for gtype in seen_gtype:
                    g = gtype.replace("-", "_")
                    if gtype in cmb_dict[table_id][s]["two"]:
                        k = "%s%s_sites"%(s_lbl,g)
                        stats_obj[table_id][k] = len(cmb_dict[table_id][s]["two"][gtype].keys())
                    if gtype in cmb_dict[table_id][s]["three"]:
                        k = "%s%s_glycan_sites"%(s_lbl,g)
                        stats_obj[table_id][k] = len(cmb_dict[table_id][s]["three"][gtype].keys())
                    if gtype in cmb_dict[table_id][s]["four"]:
                        k = "%s%s_glycans"%(s_lbl,g)
                        stats_obj[table_id][k] = len(cmb_dict[table_id][s]["four"][gtype].keys())

 
    #doc["section_stats"] = stats_obj
    doc["section_stats"] = []
    for table_id in stats_obj:
        obj = {"table_id":table_id, "table_stats":[], "sort_fields":sort_obj[table_id]}
        for f in stats_obj[table_id]:
            obj["table_stats"].append({"field":f, "count":stats_obj[table_id][f]})
        #map it to the old fields
        if table_id.find("glycosylation") != -1:
            total = total_dict[table_id] if table_id in total_dict else 0
            obj["table_stats"] = transform_section_stats(obj["table_stats"])
            obj["table_stats"] = [{"field":"total", "count":total}] + obj["table_stats"]
        doc["section_stats"].append(obj)
   
    #print (json.dumps(doc["section_stats"], indent=4))
    #exit()

 
 
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
            "biomarker_component":["total"]
            ,"citation":["total"]
        }
        ,"publication":{
            "referenced_proteins":["total"]
            ,"referenced_glycans":["total"]
            ,"glycosylation":[]
            ,"phosphorylation":[]
        }
    }
    tmp_dict = sec_map[record_type]
    sort_obj = {}
    stats_obj = {}
    for table_id in tmp_dict:
        obj = {}
        for f in tmp_dict[table_id]:
            obj[f] = 0    
        stats_obj[table_id] = obj
        sort_obj[table_id] = []


    for table_id in tmp_dict:
        if record_type == "publication" and table_id in ["glycosylation", "phosphorylation"]:
            tmp_doc = {table_id:doc[table_id]}
            get_protein_sec_stats(tmp_doc, "protein")
            stats_obj[table_id] = {}
            for obj in tmp_doc["section_stats"]:
                if obj["table_id"].split("_")[0] != table_id:
                    continue
                sort_obj[obj["table_id"]] = obj["sort_fields"]
                stats_obj[obj["table_id"]] = {}
                for o in obj["table_stats"]:
                    stats_obj[obj["table_id"]][o["field"]] = o["count"]
        elif record_type == "glycan" and table_id == "expression_tissue":
            table_id_all = table_id.split("_")[0]
            for obj in doc[table_id_all]:
                #sort_obj[table_id_all] = get_sort_fields(obj, sort_obj[table_id_all])
                get_sort_fields(obj, sort_obj[table_id_all])
                sort_obj[table_id] = sort_obj[table_id_all]
                if obj["tissue"]  != {}:
                    stats_obj[table_id]["total"] += 1
                    if "start_pos" in obj:
                        stats_obj[table_id]["total_sites"] += 1
        elif record_type == "glycan" and table_id == "expression_cell_line":
            table_id_all = table_id.split("_")[0]
            for obj in doc[table_id_all]:
                #sort_obj[table_id_all] = get_sort_fields(obj, sort_obj[table_id_all])
                get_sort_fields(obj, sort_obj[table_id_all])
                sort_obj[table_id] = sort_obj[table_id_all]
                if obj["cell_line"]  != {}: 
                    stats_obj[table_id]["total"] += 1
                    if "start_pos" in obj:
                        stats_obj[table_id]["total_sites"] += 1
        elif doc[table_id] != []:
            for obj in doc[table_id]:
                #sort_obj[table_id] = get_sort_fields(obj, sort_obj[table_id])
                get_sort_fields(obj, sort_obj[table_id])
                stats_obj[table_id]["total"] += 1
                if "start_pos" in obj:
                    stats_obj[table_id]["total_sites"] += 1
    
    #doc["section_stats"] = stats_obj
    doc["section_stats"] = []
    for table_id in stats_obj:
        if table_id == "expression_tissue":
            tmp_list = []
            for f in sort_obj[table_id]:
                if f.find("cell_line") == -1:
                    tmp_list.append(f)
            sort_obj[table_id] = tmp_list
        elif table_id == "expression_cell_line":
            tmp_list = []
            for f in sort_obj[table_id]:
                if f.find("tissue") == -1:
                    tmp_list.append(f)
            sort_obj[table_id] = tmp_list

        obj = {"table_id":table_id, "table_stats":[], "sort_fields":sort_obj[table_id]}
        if record_type == "motif" and table_id == "glycans":
            obj["table_id"] = "results"
        for f in stats_obj[table_id]:
            if f != "sort_fields":
                obj["table_stats"].append({"field":f, "count":stats_obj[table_id][f]})
        doc["section_stats"].append(obj)
    
    
    return






def update_cmb_dict(cmb_dict, table_id, site_type, start_pos, glytoucan_ac, gtype):
    if table_id not in cmb_dict:
        cmb_dict[table_id] = {
            "fuzzy":{"one":{}, "two":{}, "three":{}, "four":{}},
            "exact":{"one":{}, "two":{}, "three":{}, "four":{}}
        }
    cmb_dict[table_id][site_type]["one"][start_pos] = True
    if gtype not in cmb_dict[table_id][site_type]["two"]:
        cmb_dict[table_id][site_type]["two"][gtype] = {}
    cmb_dict[table_id][site_type]["two"][gtype][start_pos] = True
    if glytoucan_ac != "":
        if gtype not in cmb_dict[table_id][site_type]["three"]:
            cmb_dict[table_id][site_type]["three"][gtype] = {}
        cmb_dict[table_id][site_type]["three"][gtype][start_pos] = True
        cmb = "%s|%s" % (start_pos,glytoucan_ac)
        if gtype not in cmb_dict[table_id][site_type]["four"]:
            cmb_dict[table_id][site_type]["four"][gtype] = {}
        cmb_dict[table_id][site_type]["four"][gtype][cmb] = True

    return





def transform_section_stats(obj_list):

    tmp_dict = {
        "total_sites":0, 
        "n_linked_sites":0, "o_linked_sites":0, 
        "n_linked_glycan_sites":0,"o_linked_glycan_sites":0, 
        "n_linked_glycans":0, "o_linked_glycans":0,
    }
    for obj in obj_list:
        tmp_dict[obj["field"]] = obj["count"] 
        # total_sites
        if obj["field"] in ["exact_sites", "fuzzy_sites"]:
            tmp_dict["total_sites"] += obj["count"]
        #n_linked_sites
        if obj["field"] in ["exact_n_linked_sites", "fuzz_n_linked_sites"]:
            tmp_dict["n_linked_sites"] += obj["count"]
        #o_linked_sites
        if obj["field"] in ["exact_o_linked_sites", "fuzz_o_linked_sites"]:
            tmp_dict["o_linked_sites"] += obj["count"]
        #n_linked_glycan_sites
        if obj["field"] in ["exact_n_linked_glycan_sites", "fuzz_n_linked_glycan_sites"]:
            tmp_dict["n_linked_glycan_sites"] += obj["count"]
        #o_linked_glycan_sites
        if obj["field"] in ["exact_o_linked_glycan_sites", "fuzz_o_linked_glycan_sites"]:
            tmp_dict["o_linked_glycan_sites"] += obj["count"]
        #n_linked_glycans
        if obj["field"] in ["exact_n_linked_glycans", "fuzz_n_linked_glycans"]:
            tmp_dict["n_linked_glycans"] += obj["count"]
        #o_linked_glycans
        if obj["field"] in ["exact_o_linked_glycans", "fuzz_o_linked_glycans"]:
            tmp_dict["o_linked_glycans"] += obj["count"]

    tmp_list = []
    for field in tmp_dict:
        tmp_list.append({"field":field, "count":tmp_dict[field]})
    
    return tmp_list


