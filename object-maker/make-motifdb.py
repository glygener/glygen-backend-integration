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
import section_stats





##################
def get_sort_key_value_pub(obj):
    return obj["date"]

def get_sort_key_value_mut(obj):
    return obj["ann_score"]






def get_sorting_key(obj):
    return obj['sortorder']


def load_motif_glytoucan_ac_list():
    
    
    seen = {}
    data_frame = {}
    in_file = path_obj["reviewed"]+ "/glycan_motif.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        motif_ac = row[f_list.index("motif_ac")]
        motif_ac_xref = row[f_list.index("motif_ac_xref")]
        seen[motif_ac_xref] = True


    return seen.keys()









def load_property_objlist(tmp_obj_dict,in_file, prop_dict,xref_info, combo_flist_one, combo_flist_two, anchor_field):

    record_count = 0
    local_seen_dict = {}

    seen_dict = tmp_obj_dict["seen"]
    data_frame = {}
    libgly.load_sheet_as_dict(data_frame, in_file, ",", anchor_field)
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        if main_id.strip() == "":
            continue
        for tmp_row in data_frame["data"][main_id]:
            combo_id = main_id
            for f in combo_flist_one:
                combo_id += "|" + tmp_row[tmp_fl.index(f)]
            combo_id = combo_id.strip()
            if combo_id not in local_seen_dict:
                record_count += 1
                local_seen_dict[combo_id] = True

            obj_one = {}    
            for prop in prop_dict:
                f = prop_dict[prop]
                obj_one[prop] = tmp_row[tmp_fl.index(f)]
                if prop == "do_id":
                    xref_key = tmp_row[tmp_fl.index(xref_info[0])]
                    xref_id = tmp_row[tmp_fl.index(xref_info[1])]
                    do_id = combo_id.split("|")[1]
                    obj_one[prop] = do_id
                    if do_id == "":
                        combo_id = "%s|%s-%s" % (main_id, xref_key.split("_")[-1], xref_id)
                        obj_one[prop] = "%s-%s" % (xref_key.split("_")[-1], xref_id)
                        database_label = tmp_row[tmp_fl.index("database_label")]
                        do_name = "%s [%s disease name]" % (database_label, xref_badge)
                        obj_one["name"] = do_name
                        obj_one["url"] = libgly.get_xref_url(map_dict, "protein_xref_do_placeholder", "",is_cited)
                    else:
                        do_name = doid2name[do_id][0]
                        obj_one["name"] = do_name[0].upper() + do_name[1:] + " [DO disease name]"
                        obj_one["url"] = libgly.get_xref_url(map_dict, "protein_xref_do", do_id,is_cited)
                        if "protein_xref_icd10cm" in doid2xrefid[do_id]:
                            obj_one["icd10"] = doid2xrefid[do_id]["protein_xref_icd10cm"][0]

            if combo_id not in seen_dict:
                seen_dict[combo_id] = True
                tmp_obj_dict[combo_id] = obj_one
                if combo_flist_two != []:
                    obj_one["evidence"] = []


            if combo_flist_two != []:
                xref_key = tmp_row[tmp_fl.index(xref_info[0])]
                xref_id = tmp_row[tmp_fl.index(xref_info[1])]
                xref_id = main_id if xref_id.strip() == "" else xref_id
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                obj_two = {"database":xref_badge, "id":xref_id, "url":xref_url}
                combo_id_xref = combo_id
                for f in combo_flist_two:
                    combo_id_xref += "|" + tmp_row[tmp_fl.index(f)]
                combo_id_xref = combo_id_xref.strip()
                if combo_id_xref not in seen_dict:
                    seen_dict[combo_id_xref] = True
                    tmp_obj_dict[combo_id]["evidence"].append(obj_two)
    
    return record_count



def load_properity_list(tmp_obj_dict, in_file, field_list, sep):

    data_frame = {}
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        if main_id not in tmp_obj_dict:
            tmp_obj_dict[main_id] = []
        for tmp_row in data_frame["data"][main_id]:
            for k in field_list:
                field_value = tmp_row[tmp_fl.index(k)].strip()
                value_list = [field_value] if sep == "" else field_value.split(sep)
                for value in value_list:
                    if value != "" and value not in tmp_obj_dict[main_id]:
                        tmp_obj_dict[main_id].append(value.strip())

    return






def update_record_stat(record_stat, file_name, prop_name, n, combo_list):

    if file_name not in record_stat:
        record_stat[file_name] = {"recordfields":combo_list}

    if prop_name not in record_stat[file_name]:
        record_stat[file_name][prop_name] = 0

    record_stat[file_name][prop_name] += n
    return


#######################################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-s","--sec",action="store",dest="sec",help="Object section")

    (options,args) = parser.parse_args()
    sec_name = options.sec


    global config_obj
    global path_obj
    global map_dict
    global data_dir
    global main_dict


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"


    global is_cited
    is_cited = libgly.get_is_cited()


    expected_dslist = []
    sec_info = json.loads(open("generated/misc/motif_sectioninfo.json", "r").read())
    expected_dslist =  csvutil.get_expected_dslist("motif", [], {})
    #print (json.dumps(expected_dslist, indent=4))
    #print (len(expected_dslist))
    #exit()

    residue_heirarchy = json.loads(open("generated/misc/residue_heirarchy.json", "r").read())
    residue2class = {}
    for c in residue_heirarchy:
        for m in residue_heirarchy[c]:
            residue2class[m] = c
    
    main_dict = {}
    for sec in sec_info:
        main_dict[sec] = {"seen":{}}

    map_dict = {}
    csvutil.load_dictionaries(map_dict, "generated/misc/")
   

    motif_glytoucan_ac_list = load_motif_glytoucan_ac_list()
    record_stat = {}
    glycan_type_dict = {}
    motif_name_dict = {}
    main_obj_dict = {}
   


    biomarker_dict = csvutil.get_biomarker_dict("glycan") 
    

    log_file = "logs/make-motifdb.log"
    msg = "make-motifdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")

    file_idx, file_count = 1, len(expected_dslist) 
    for file_name in expected_dslist:
        in_file = "%s%s.csv" % (data_dir,file_name)
        if os.path.isfile(in_file) == False:
            msg = "make-motifdb: file %s does NOT exist!" % (in_file)
            csvutil.write_log_msg(log_file, msg, "a")
            sys.exit()

        msg = "make-motifdb: %s [%s/%s]" % (in_file, file_idx, file_count)
        csvutil.write_log_msg(log_file, msg, "a")
        file_idx += 1

        #--> glytoucan_ac
        sheet_name = "glycan_masterlist"
        if file_name.find(sheet_name) != -1:
            species = file_name.split("_")[0]
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
            tmp_fl = data_frame["fields"]
            for main_id in data_frame["data"]:
                if main_id in main_dict["glytoucan_ac"]:
                    continue
                main_dict["glytoucan_ac"][main_id] = main_id
                tmp_row = data_frame["data"][main_id][0]
                mass = tmp_row[tmp_fl.index("glycan_mass")]
                mass_pme = tmp_row[tmp_fl.index("glycan_permass")]
                number_monosaccharides = tmp_row[tmp_fl.index("monosaccharides")]

                glycan_type_dict[main_id] = tmp_row[tmp_fl.index("glytoucan_type")]
                if mass != "":
                    prop_name = "mass"
                    mass = round(float(tmp_row[tmp_fl.index("glycan_mass")]), 2)
                    main_dict[prop_name][main_id] = mass
                    combo_list = ["glytoucan_ac", "mass"]
                    update_record_stat(record_stat, file_name, prop_name, 1, combo_list)
                if mass_pme != "":
                    prop_name = "mass_pme"
                    mass_pme = round(float(tmp_row[tmp_fl.index("glycan_permass")]), 2)
                    main_dict[prop_name][main_id] = mass_pme
                    combo_list = ["glytoucan_ac", "mass_pme"]
                    update_record_stat(record_stat, file_name, prop_name, 1,combo_list)
                if number_monosaccharides.isdigit():
                    prop_name = "number_monosaccharides"
                    number_monosaccharides = int(tmp_row[tmp_fl.index("monosaccharides")])
                    main_dict[prop_name][main_id] = number_monosaccharides
                    combo_list = ["glytoucan_ac", "number_monosaccharides"]
                    update_record_stat(record_stat, file_name, prop_name, 1,combo_list)
        
        #--> fully_determined
        sheet_name = "fully_determined"
        if file_name.find(sheet_name) != -1:
            species = file_name.split("_")[0]
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
            tmp_fl = data_frame["fields"]
            for main_id in data_frame["data"]:
                main_dict["fully_determined"][main_id] = "yes"

        #--> sequences
        sheet_name = "_sequences_"
        if file_name.find(sheet_name) != -1:
            species = file_name.split("_")[0]
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
            tmp_fl = data_frame["fields"]
            prop2field = {
                "iupac":"sequence_iupac_extended",
                "wurcs":"sequence_wurcs",
                "glycoct":"sequence_glycoct",
                "inchi":"sequence_inchi",
                "smiles_isomeric":"sequence_smiles_isomeric",
                "glycam":"sequence_glycam_iupac"
            }
            for main_id in data_frame["data"]:
                tmp_row = data_frame["data"][main_id][0]
                for p in prop2field:
                    field = prop2field[p]
                    if field in tmp_fl:
                        main_dict[p][main_id] = tmp_row[tmp_fl.index(field)]
                        combo_list = ["glytoucan_ac", field]
                        update_record_stat(record_stat, file_name, p, 1,combo_list)
                        if field == "sequence_inchi":
                            inchi_key = tmp_row[tmp_fl.index("inchi_key")]
                            url = libgly.get_xref_url(map_dict, "glycan_xref_inchi_key", inchi_key, is_cited)
                            o = {"key":inchi_key, "url":url}
                            main_dict["inchi_key"][main_id] = o
                            combo_list = ["glytoucan_ac", field]
                            update_record_stat(record_stat, file_name, "inchi_key", 1, combo_list)
        
        #--> biomarkers
        for sheet_name in ["glycan_biomarkers"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "biomarkers"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    #print (main_id, main_id in biomarker_dict)
                    if main_id in biomarker_dict:
                        for doc in biomarker_dict[main_id]:
                            combo_id = "%s|%s" % (main_id, doc["biomarker_id"])
                            if "aclist" in doc:
                                doc.pop("aclist")
                            load_obj[combo_id] = doc
 
        
        #--> dictionary
        for sheet_name in ["glycan_dictionary"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "dictionary"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        o = {}
                        for f in tmp_fl:
                            if f in ["xref_id", "xref_key"]:
                                continue
                            o[f] = tmp_row[tmp_fl.index(f)]
                        o["synonymns"] = o["synonymns"].split("|")  
                        xref_id = tmp_row[tmp_fl.index("xref_id")]
                        xref_key = tmp_row[tmp_fl.index("xref_key")]
                        xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                        o["evidence"] = [{"id":xref_id, "url":xref_url, "database":xref_badge}]
                        combo_id = "%s|%s|%s" % (main_id, xref_id, o["term"])
                        load_obj[combo_id] = o



        #--> publication
        for sheet_name in ["glycan_citations_motif"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "publication"
                prop_dict = {"title":"title", 
                        "journal":"journal_name","date":"publication_date","authors":"authors"}
                xref_info = ["src_xref_key", "src_xref_id"]
                combo_flist_one = ["xref_key", "xref_id"]
                combo_flist_two = ["src_xref_key", "src_xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj,in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two, "glytoucan_ac") 
                combo_list = ["glytoucan_ac"] + combo_flist_one
                update_record_stat(record_stat, file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id, xref_key, xref_id = combo_id.split("|")
                    xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    load_obj[combo_id]["date"] = load_obj[combo_id]["date"].split(" ")[0]
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    load_obj[combo_id]["reference"] = [
                        {
                            "type":xref_badge,
                            "id":xref_id,
                            "url":xref_url
                        }
                    ]
                    


        #--> glycans
        for sheet_name in ["glycan_motif"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "glycans"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        motif_ac = tmp_row[tmp_fl.index("motif_ac")]
                        motif_ac_xref = tmp_row[tmp_fl.index("motif_ac_xref")]
                        motif_name = tmp_row[tmp_fl.index("motif_name")]
                        combo_id = "%s|%s" % (motif_ac_xref,main_id)
                        xref_key = "glycan_xref_glytoucan"
                        xref_url = libgly.get_xref_url(map_dict, xref_key, main_id,is_cited)
                        main_dict["glycans"][combo_id] = {"glytoucan_ac":main_id, 
                                "url":xref_url
                        }
                        if motif_ac_xref not in main_dict["motif_ac"]:
                            main_dict["motif_ac"][motif_ac_xref] = []
                        if motif_ac not in main_dict["motif_ac"][motif_ac_xref]:
                            main_dict["motif_ac"][motif_ac_xref].append(motif_ac)
                        xref_url = libgly.get_xref_url(map_dict, xref_key, motif_ac_xref,is_cited)
                        xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                        o = {"id":motif_ac_xref, "database": xref_badge, "url":xref_url }
                        combo_id = "%s|%s|%s" % (motif_ac_xref,xref_key, motif_ac_xref)
                        main_dict["crossref"][combo_id] = o

                        name, domain = motif_name, "motifname"
                        #combo_id = "%s|%s|%s" % (motif_ac_xref, name, domain)
                        #main_dict["names"][combo_id] = {"name":name, "domain":domain}
                        if motif_ac not in motif_name_dict:
                            motif_name_dict[motif_ac] = {}
                        motif_name_dict[motif_ac][name] = True
        #--> names
        for sheet_name in ["glycan_names"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "names"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        name = tmp_row[tmp_fl.index("glycan_name")]
                        domain = tmp_row[tmp_fl.index("glycan_name_domain")]
                        combo_id = "%s|%s|%s" % (main_id, name, domain)
                        load_obj[combo_id] = {"name":name, "domain":domain}

        
        #--> synonym
        for sheet_name in ["glycan_motif"]:
            if file_name.find(sheet_name) != -1:
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "motif_ac_xref")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    seen_kw = {}
                    for tmp_row in data_frame["data"][main_id]:
                        val_dict = {
                            "synonym":tmp_row[tmp_fl.index("alternative_name")],
                            "keywords": tmp_row[tmp_fl.index("keyword")],
                            "reducing_end": tmp_row[tmp_fl.index("reducing_end")],
                            "aglycon":  tmp_row[tmp_fl.index("aglycon")],
                            "alignment_method": tmp_row[tmp_fl.index("alignment")]
                        }
                        for p in val_dict:
                            val = val_dict[p].strip()
                            if val == "":
                                continue
                            for val in val.split(";"):
                                val = val.strip()
                                if p not in main_dict:
                                    main_dict[p] = {}
                                if p in ["reducing_end", "aglycon", "alignment_method"]:
                                    main_dict[p][main_id] = val
                                elif p in ["synonym", "keywords"]:
                                    if main_id not in main_dict[p]:
                                        main_dict[p][main_id] = []
                                    if val not in main_dict[p][main_id]:
                                        if p == "keywords":
                                            xref_key = "glycan_xref_motif"
                                            xref_id = val.replace(" ", "_")
                                            xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                                            keyword_obj = {"label":val, "url":xref_url}
                                            if val not in seen_kw:
                                                main_dict[p][main_id].append(keyword_obj)
                                                seen_kw[val] = True
                                        else:
                                            main_dict[p][main_id].append(val)
                            

        #--> subsumption
        for sheet_name in ["glycan_subsumption"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "subsumption"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        related_accession = tmp_row[tmp_fl.index("related_accession")]
                        if related_accession not in glycan_type_dict:
                            continue
                        glytoucan_type = glycan_type_dict[related_accession]
                        combo_id = "%s|%s" % (main_id, related_accession)
                        load_obj[combo_id] = {"id":related_accession, "type":glytoucan_type}



        #--> residues
        for sheet_name in ["glycan_enzyme"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "residues"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    r_count = 1
                    for tmp_row in data_frame["data"][main_id]:
                        canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                        residue_id = tmp_row[tmp_fl.index("residue_id")]
                        residue_name = tmp_row[tmp_fl.index("residue_name")]
                        parent_residue_id = tmp_row[tmp_fl.index("parent_residue_id")]
                        attached_by = "rxn.%s" % (canon)
                        detached_by = ""
                        r = {"id":residue_id, "name":residue_name, 
                                "attachedby":attached_by, 
                                "detachedby":detached_by, "parentid":parent_residue_id}
                        combo_id = "%s|%s" % (main_id, r_count)
                        load_obj[combo_id] = r
                        r_count += 1
                combo_list = ["glytoucan_ac", "residue_count"] 
                update_record_stat(record_stat, file_name, prop_name, len(load_obj.keys()), combo_list)



        #--> crossref
        #--> In addition to GlyToucan crossref, also do dictionary xref
        for sheet_name in ["glycan_xref_dictionary"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "crossref"
                load_obj = main_dict[prop_name]
                prop_dict = {"id":"xref_id"}
                combo_flist_one = ["xref_key", "xref_id"]
                xref_info = []
                combo_flist_two = []
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two, "glytoucan_ac")
                combo_list = ["glytoucan_ac"] + combo_flist_one
                update_record_stat(record_stat, file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id = combo_id.split("|")[0]
                    uniprotkb_ac = main_id.split("-")[0]
                    xref_key, xref_id = combo_id.split("|")[-2], combo_id.split("|")[-1]
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    load_obj[combo_id]["url"] = xref_url
                    load_obj[combo_id]["database"] = xref_badge
                    cats_dict = map_dict["xrefkey2category"]
                    xref_categories = cats_dict[xref_key] if xref_key in cats_dict else []
                    load_obj[combo_id]["categories"] = xref_categories

 
    
        #--> composition (both composition and composition_advanced used advanced dataset)              
        for sheet_name in ["glycan_monosaccharide_composition_advanced"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "composition_expanded"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    type_combo_list = []
                    for tmp_row in data_frame["data"][main_id]:
                        for residue in map_dict["residue2name"]:
                            if residue not in tmp_fl:
                                continue
                            n = int(tmp_row[tmp_fl.index(residue)])
                            name = map_dict["residue2name"][residue][0]
                            cid = map_dict["residue2cid"][residue][0]
                            residue = "other" if residue.lower() == "xxx" else residue.lower()
                            o = {"name":name, "residue":residue, "count":n}
                            if cid != "" and cid.strip() != "0":
                                url = libgly.get_xref_url(map_dict, "glycan_xref_monosaccharide_residue_name", cid,is_cited)
                                o = {"name":name, "residue":residue, "count":n, "cid":cid, "url":url}
                            combo_id = "%s|%s|%s" % (main_id, residue, cid)
                            load_obj[combo_id] = o
                combo_list = ["glytoucan_ac", "residue"]
                update_record_stat(record_stat, file_name, prop_name, len(load_obj.keys()), combo_list)

        #--> composition_advanced               
        for sheet_name in ["glycan_monosaccharide_composition_advanced"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "composition"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    type_combo_list = []
                    for tmp_row in data_frame["data"][main_id]:
                        class_count = {}
                        for residue in map_dict["residue2name"]:
                            if residue not in tmp_fl:
                                continue
                            #consider only top level residues
                            if residue not in residue_heirarchy:
                                continue
                            r_class = residue
                            n = int(tmp_row[tmp_fl.index(residue)])
                            class_count[r_class] = n
                        
                        for r_class in class_count:
                            n = class_count[r_class]
                            r_class_name = map_dict["residue2name"][r_class][0]
                            cid = map_dict["residue2cid"][r_class][0]
                            r_class = "other" if r_class.lower() == "xxx" else r_class.lower()
                            o = {"name":r_class_name, "residue":r_class, "count":n}
                            if cid != "" and cid.strip() != "0":
                                url = libgly.get_xref_url(map_dict, "glycan_xref_monosaccharide_residue_name",cid,is_cited)
                                o = {"name":r_class_name, "residue":r_class, "count":n, "cid":cid, "url":url}
                            combo_id = "%s|%s|%s" % (main_id, r_class, cid)
                            load_obj[combo_id] = o
                combo_list = ["glytoucan_ac", "residue"]
                update_record_stat(record_stat, file_name, prop_name, len(load_obj.keys()), combo_list)









    tmp_dict = {}
    for main_id in main_dict["glytoucan_ac"]:
        tmp_dict[main_id] = {}
        for sec in sec_info:
            if sec_info[sec]["category"] in ["string"]:
                tmp_dict[main_id][sec] = ""
            elif sec_info[sec]["category"] in ["float"]:
                tmp_dict[main_id][sec] = 0.0
            elif sec_info[sec]["category"] in ["int"]:
                tmp_dict[main_id][sec] = 0
            elif sec_info[sec]["category"] in ["obj"]:
                tmp_dict[main_id][sec] = {}
            elif sec_info[sec]["category"] in ["list", "objlist"]:
                tmp_dict[main_id][sec] = []



    for sec in main_dict:
        if type(main_dict[sec]) is dict:
            if "seen" in main_dict[sec]:
                main_dict[sec].pop("seen")
        for combo_id in main_dict[sec]:
            main_id = combo_id.split("|")[0]
            if main_id not in tmp_dict:
                continue
            if sec_info[sec]["category"] in ["string", "int", "float", "list", "obj"]:
                tmp_dict[main_id][sec] =  main_dict[sec][combo_id]
            elif sec_info[sec]["category"] in ["objlist"]:
                o = main_dict[sec][combo_id]
                if "deleteflag" not in o:
                    tmp_dict[main_id][sec].append(o)
        


    record_count = 0
    for main_id in tmp_dict:
        if main_id == "seen":
            continue
        obj = tmp_dict[main_id]
        if main_id not in motif_glytoucan_ac_list:
            continue
        if "motif_ac" not in obj:
            continue
        motif_ac_list = []
        for motif_ac in obj["motif_ac"]:
            motif_ac_list.append(motif_ac)
        obj.pop("motif_ac")
        parent_names = obj["names"]
        for motif_ac in motif_ac_list:
            obj["glytoucan_ac"] = main_id
            obj["motif_ac"] = motif_ac
            motif_names = []
            if motif_ac in motif_name_dict:
                for motif_name in motif_name_dict[motif_ac]:
                    motif_names.append({"name":motif_name, "domain":"motifname"})
            obj["names"] = parent_names + motif_names
            
            section_stats.get_sec_stats(obj, "motif")
            out_file = "jsondb/motifdb/%s.json" % (motif_ac)
            for k in ["mass", "mass_pme", "number_monosaccharides"]:
                if obj[k] == 0.0 or obj[k] == 0:
                    obj.pop(k)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(obj, indent=4)))
            record_count += 1 
    msg = "make-motifdb: final filtered in: %s motif objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "a")






if __name__ == '__main__':
        main()



