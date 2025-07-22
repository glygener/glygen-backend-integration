#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import csvutil
import libgly


def load_disease2glycan_dict():

    seen = {}
    file_list = glob.glob("reviewed/glycan_biomarkers.csv")

    for in_file in file_list:
        if in_file.find(".stat.csv") != -1:
            continue
        data_frame = {}
        csvutil.load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            gtc = row[f_list.index("glytoucan_ac")]
            do_id = row[f_list.index("do_id")] if "do_id" in f_list else ""
            mondo_id = row[f_list.index("mondo_id")] if "mondo_id" in f_list else ""
            mim_id = row[f_list.index("mim_id")] if "mim_id" in f_list else ""
            disease_id = "doid.%s" % (do_id) if do_id != "" else ""
            xref_id = row[f_list.index("xref_id")] if "xref_id" in f_list else ""
            xref_key = row[f_list.index("xref_key")] if "xref_key" in f_list else ""
            
            if disease_id == "":
                disease_id = "mondo.%s" % (mondo_id) if mondo_id != "" else ""
            if disease_id == "":
                disease_id = "mim.%s" % (mim_id) if mim_id != "" else ""
            if disease_id == "":
                continue
            if disease_id not in seen:
                seen[disease_id] = {}
            if gtc not in seen[disease_id]:
                seen[disease_id][gtc] = {}
            if xref_key not in seen[disease_id][gtc]:
                seen[disease_id][gtc][xref_key] = {}
            seen[disease_id][gtc][xref_key][xref_id] = True


    return seen


def load_disease2biomarker_dict():

    seen = {}
    file_list = glob.glob("reviewed/glycan_biomarkers.csv")
    file_list += glob.glob("reviewed/*_protein_biomarkers.csv")

    for in_file in file_list:
        if in_file.find(".stat.csv") != -1:
            continue
        data_frame = {}
        csvutil.load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            biomarker_id = row[f_list.index("biomarker_id")]
            do_id = row[f_list.index("do_id")] if "do_id" in f_list else ""
            mondo_id = row[f_list.index("mondo_id")] if "mondo_id" in f_list else ""
            mim_id = row[f_list.index("mim_id")] if "mim_id" in f_list else ""
            disease_id = "doid.%s" % (do_id) if do_id != "" else ""
            xref_id = row[f_list.index("xref_id")] if "xref_id" in f_list else ""
            xref_key = row[f_list.index("xref_key")] if "xref_key" in f_list else ""
         
            if disease_id == "":
                disease_id = "mondo.%s" % (mondo_id) if mondo_id != "" else ""
            if disease_id == "":
                disease_id = "mim.%s" % (mim_id) if mim_id != "" else ""
            if disease_id == "":
                continue
            if disease_id not in seen:
                seen[disease_id] = {}
            if biomarker_id not in seen[disease_id]:
                seen[disease_id][biomarker_id] = {}
            if xref_key not in seen[disease_id][biomarker_id]:
                seen[disease_id][biomarker_id][xref_key] = {}
            seen[disease_id][biomarker_id][xref_key][xref_id] = True

    return seen


def load_disease2protein_dict():

    seen = {}
    file_list = glob.glob("reviewed/*_protein_disease_*.csv")
    file_list += glob.glob("reviewed/*_protein_biomarkers.csv")
    #file_list += glob.glob("reviewed/*_protein_mutation_cancer.csv")
    #file_list += glob.glob("reviewed/*_protein_mutation_somatic.csv")
    #file_list += glob.glob("reviewed/*_protein_mutation_germline.csv")
    #file_list += glob.glob("reviewed/*_protein_expression_disease.csv")
    
    for in_file in file_list:
        if in_file.find(".stat.csv") != -1:
            continue
        data_frame = {}
        csvutil.load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            do_id = row[f_list.index("do_id")] if "do_id" in f_list else ""
            xref_id = row[f_list.index("xref_id")] if "xref_id" in f_list else ""
            xref_key = row[f_list.index("xref_key")] if "xref_key" in f_list else ""
            mondo_id = row[f_list.index("mondo_id")] if "mondo_id" in f_list else ""
            mim_id = row[f_list.index("mim_id")] if "mim_id" in f_list else ""
            disease_id = "doid.%s" % (do_id) if do_id != "" else ""
            if disease_id == "":
                disease_id = "mondo.%s" % (mondo_id) if mondo_id != "" else ""
            if disease_id == "":
                disease_id = "mim.%s" % (mim_id) if mim_id != "" else ""
            if disease_id == "":
                continue
            if disease_id not in seen:
                seen[disease_id] = {}
            if canon not in seen[disease_id]:
                seen[disease_id][canon] = {}
            if xref_key not in seen[disease_id][canon]:
                seen[disease_id][canon][xref_key] = {}
            seen[disease_id][canon][xref_key][xref_id] = True

    return seen


def get_glycan_objects(in_dict):

    obj_list = []
    idx = 0
    total = len(in_dict.keys())
    for gtc in in_dict:
        idx += 1
        in_file = "jsondb/glycandb/%s.json" % (gtc)
        doc = json.load(open(in_file))
        obj = {}
        for p in ["glytoucan_ac", "names"]:
            obj[p] = doc[p]

        obj["species"] = []
        for o in doc["species"]:
            oo = {}
            for p in ["name", "taxid", "common_name", "glygen_name", "reference_species"]:
                oo[p] = o[p]
            obj["species"].append(oo)
        
        obj["evidence"] = []
        for xref_key in in_dict[gtc]:
            for xref_id in in_dict[gtc][xref_key]:
                xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                o = {"id":xref_id, "database":xref_badge, "url":xref_url}
                obj["evidence"].append(o)
        obj_list.append(obj)

    return obj_list


def get_biomarker_objects(in_dict):

    p_list_one = ["biomarker_id", "biomarker_canonical_id", "best_biomarker_role"]
    p_list_two = ["assessed_biomarker_entity", "biomarker", "assessed_biomarker_entity_id"]
    p_list_two += ["assessed_entity_type"]

    
    obj_list = []
    idx = 0
    total = len(in_dict.keys())
    for biomarker_id in in_dict:
        idx += 1
        in_file = "jsondb/biomarkerdb/%s.json" % (biomarker_id)
        doc = json.load(open(in_file))
        obj = {}
        for p in p_list_one:
            obj[p] = doc[p]
        #obj["biomarker_component"] = []
        #for cmp_obj in doc["biomarker_component"]:
        #    o = {}
        #    for p in p_list_two:
        #        o[p] = cmp_obj[p]
        #    o["evidence"] = []
        #    for oo in cmp_obj["evidence_source"]:
        #        for kk in ["evidence_list", "tags"]:
        #            if kk in oo:
        #                oo.pop(kk)
        #        o["evidence"].append(oo)
        #    obj["biomarker_component"].append(o)
        
        obj["evidence"] = []
        for xref_key in in_dict[biomarker_id]:
            for xref_id in in_dict[biomarker_id][xref_key]:
                xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                o = {"id":xref_id, "database":xref_badge, "url":xref_url}
                obj["evidence"].append(o)
        obj_list.append(obj)

    return obj_list




def get_protein_objects(in_dict):

    obj_list = []
    idx = 0
    total = len(in_dict.keys())
    for canon in in_dict:
        idx += 1
        in_file = "jsondb/proteindb/%s.json" % (canon)
        doc = json.load(open(in_file))
        obj = {}
        for p in ["uniprot_canonical_ac", "uniprot_ac", "uniprot_id"]:
            obj[p] = doc[p]
        obj["refseq_ac"] = ""
        if "refseq" in doc:
            if "ac" in doc["refseq"]:
                obj["refseq_ac"] = doc["refseq"]["ac"]
        obj["glycoprotein"] = "no"
        if "keywords" in doc:
            if "glycoprotein" in doc["keywords"]:
                obj["glycoprotein"] = "yes"
        
        obj["species"] = []
        for o in doc["species"]:
            oo = {}
            for p in ["name", "taxid", "common_name", "glygen_name", "reference_species"]:
                oo[p] = o[p]
            obj["species"].append(oo) 
        for k in ["protein_names", "gene_names"]:
            obj[k] = []
            for o in doc[k]:
                oo = {}
                for p in ["name", "resource", "type"]:
                    oo[p] = o[p]
                obj[k].append(oo)
        obj["evidence"] = []
        for xref_key in in_dict[canon]:
            for xref_id in in_dict[canon][xref_key]:
                xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                o = {"id":xref_id, "database":xref_badge, "url":xref_url}
                obj["evidence"].append(o)
        obj_list.append(obj) 


    return obj_list


def get_xref_obj_list(doc):

    seen = {}
    obj_list = []
    obj = doc["recommended_name"]
    xref_id, xref_db, xref_url = obj["id"], obj["resource"], obj["url"]
    obj_list.append({"id":xref_id, "database":xref_db, "url":xref_url, "categories":[]})
    seen[xref_id] = True
    for obj in doc["synonyms"]:
        xref_id, xref_db, xref_url = obj["id"], obj["resource"], obj["url"]
        if xref_id not in seen:
            obj_list.append({"id":xref_id, "database":xref_db, "url":xref_url, "categories":[]})
            seen[xref_id] = True

    return obj_list




def main():

    global config_obj
    global path_obj
    global map_dict
    global is_cited




    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]


    is_cited = libgly.get_is_cited()




    log_file = "logs/update-diseasedb.log"
    msg = "update-diseasedb: started"
    csvutil.write_log_msg(log_file, msg, "w")

    map_dict = {}
    csvutil.load_dictionaries(map_dict, "generated/misc/")



    file_list = glob.glob("jsondb/interm/diseasedb/*.json")
    #file_list = glob.glob("jsondb/interm/diseasedb/doid.10534.json")
    #file_list += glob.glob("jsondb/interm/diseasedb/doid.1612.json")

    disease2protein = load_disease2protein_dict()
    disease2glycan = load_disease2glycan_dict()
    disease2biomarker = load_disease2biomarker_dict()


    record_count = 0
    for in_file in file_list:
        doc = json.load(open(in_file))
        disease_id = in_file.split("/")[-1].replace(".json", "")
        doc["crossref"] = get_xref_obj_list(doc)
        doc["proteins"], doc["glycans"], doc["biomarkers"] = [], [], []
        if disease_id in disease2protein:
            doc["proteins"] = get_protein_objects(disease2protein[disease_id])
        if disease_id in disease2glycan:
            doc["glycans"] = get_glycan_objects(disease2glycan[disease_id])
        if disease_id in disease2biomarker:
            doc["biomarkers"] = get_biomarker_objects(disease2biomarker[disease_id])
        out_file = "jsondb/diseasedb/%s.json" % (disease_id)
        doc["record_id"] = disease_id
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        record_count += 1
        if record_count%100 == 0:
            msg = "update-diseasedb: updated %s out of %s disease records" % (record_count, len(file_list))
            csvutil.write_log_msg(log_file, msg, "a")
        
    log_file = "logs/update-diseasedb.log"
    msg = "update-diseasedb: final updated: %s disease objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "a")



if __name__ == '__main__':
    main()


