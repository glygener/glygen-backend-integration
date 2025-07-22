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
#import commands
import subprocess

import libgly
import csvutil
import batchutil
import section_stats






##################
def get_sort_key_value_pub(obj):
    return obj["date"]

def get_sort_key_value_mut(obj):
    return obj["ann_score"]






def get_sorting_key(obj):
    return obj['sortorder']


def load_subsumption_heirarchy(in_file):

    heirarchy_dict = {}
    data_frame = {}
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        for tmp_row in data_frame["data"][main_id]:
            related_accession = tmp_row[tmp_fl.index("related_accession")]
            relationship = tmp_row[tmp_fl.index("relationship")].lower()
            glytoucan_type = tmp_row[tmp_fl.index("glytoucan_type")].lower()
            gt_list = ["topology", "composition", "basecomposition"]
            rl_list = ["ancestor", "descendant"]
            if relationship in rl_list:
                if main_id not in heirarchy_dict:
                    heirarchy_dict[main_id] = {}
                if related_accession not in heirarchy_dict[main_id]:
                    heirarchy_dict[main_id][related_accession] = {}
                heirarchy_dict[main_id][related_accession] = relationship

    return heirarchy_dict

    

def load_motif2gtc():
    
    seen = {}
    data_frame = {}
    in_file = path_obj["reviewed"]+ "/glycan_motif.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")]
        motif_ac_xref = row[f_list.index("motif_ac_xref")]
        if motif_ac_xref not in seen:
            seen[motif_ac_xref] = {}
        seen[motif_ac_xref][glytoucan_ac] = True


    return seen


def load_species_info(species_obj, species_list):

    seen = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)
    species_obj_list = []
    for k in sorted(species_obj, reverse=True):
        obj = species_obj[k]
        if obj["short_name"] not in seen and obj["is_reference"] == "yes":
            o = {"shortname":obj["short_name"], "sortorder":obj["sort_order"]}
            species_obj_list.append(o)
            seen[obj["short_name"]] = True

    species_obj_list.sort(key=get_sorting_key)
    for o in species_obj_list:
        species_list.append(o["shortname"])


    return








def load_property_objlist(tmp_obj_dict,in_file, prop_dict,xref_info, combo_flist_one, combo_flist_two, anchor_field):

    record_count = 0
    local_seen_dict = {}

    seen_dict = tmp_obj_dict["seen"]


    FR = open(in_file, "r")
    idx = 0
    tmp_fl = []
    for line in FR:
        idx += 1
        tmp_row = line.strip().split("\",\"")
        tmp_row[0], tmp_row[-1] = tmp_row[0].replace("\"", ""), tmp_row[-1].replace("\"", "")
        if idx == 1:
            tmp_fl = tmp_row
            continue
        else:
            main_id = tmp_row[tmp_fl.index(anchor_field)]    
            if main_id.strip() == "":
                continue
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
                rr_list = [[xref_info[0], xref_info[1]]]
                if len(xref_info) == 4 :
                    rr_list.append([xref_info[2], xref_info[3]])
                for rr in rr_list:
                    if rr[0] not in tmp_fl or rr[1] not in tmp_fl:
                        continue
                    xref_key = tmp_row[tmp_fl.index(rr[0])]
                    xref_id = tmp_row[tmp_fl.index(rr[1])]
                    xref_id = main_id if xref_id.strip() == "" else xref_id
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    obj_two = {"database":xref_badge, "id":xref_id, "url":xref_url}
                    combo_id_xref = combo_id + "|" + xref_key + "|" + xref_id
                    #for f in combo_flist_two:
                    #    combo_id_xref += "|" + tmp_row[tmp_fl.index(f)]
                    #combo_id_xref = combo_id_xref.strip()
                    if combo_id_xref not in seen_dict:
                        seen_dict[combo_id_xref] = True
                        tmp_obj_dict[combo_id]["evidence"].append(obj_two)
    FR.close()

 
    return record_count

def load_motif_names():

    tmp_dict = {}
    in_file = "reviewed/glycan_motif.csv"
    data_frame = {}
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        for tmp_row in data_frame["data"][main_id]:
            motif_ac = tmp_row[tmp_fl.index("motif_ac")] 
            motif_ac_xref = tmp_row[tmp_fl.index("motif_ac_xref")]
            motif_name = tmp_row[tmp_fl.index("motif_name")]
            motif_syn = tmp_row[tmp_fl.index("alternative_name")]
            if main_id == motif_ac_xref:
                tmp_dict[motif_ac_xref] = {"name":motif_name, "synonym_list":motif_syn.split(";")}

    return tmp_dict



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

def get_protein_info(protein_obj, info_type):

    ret_val = ""
    if info_type == "gene_name":
        if protein_obj["gene"] != []:
            ret_val = protein_obj["gene"][0]["name"]
    elif info_type == "gene_url":
        if protein_obj["gene"] != []:
            ret_val = protein_obj["gene"][0]["url"]
    elif info_type == "tax_id":
        if protein_obj["species"] != []:
            ret_val = protein_obj["species"][0]["taxid"]
    elif info_type == "tax_name":
        if protein_obj["species"] != []:
            ret_val = protein_obj["species"][0]["name"]
    elif info_type == "tax_common_name":
        if protein_obj["species"] != []:
            if "common_name" in protein_obj["species"][0]:
                ret_val = protein_obj["species"][0]["common_name"]
    elif info_type == "recommendedname":
        if "protein_names" in protein_obj:
            syn_list = []
            for o in protein_obj["protein_names"]:
                if o["type"] == "recommended":
                    ret_val = o["name"]
                else:
                    syn_list.append(o["name"])
            if ret_val == "" and syn_list != []:
                ret_val = syn_list[0]
                
    return ret_val



def load_protein_objects(protein_obj_dict):

    file_list = glob.glob(data_dir + "/*_protein_masterlist.csv")
    for in_file in file_list:
        species = in_file.split("/")[-1].split("_")[0]
        data_frame = {}
        libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
        tmp_fl = data_frame["fields"]
        for canon in data_frame["data"]:
            if canon not in protein_obj_dict:
                protein_obj_dict[canon] = {"gene":[], "species":[]}
            uniprotkb_ac = canon.split("-")[0]
            tmp_row = data_frame["data"][canon][0]
            gene_name = tmp_row[tmp_fl.index("gene_name")]
            xref_key, xref_id = "protein_xref_uniprotkb", uniprotkb_ac
            gene_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
            o = { "name":gene_name ,"url":gene_url}
            protein_obj_dict[canon]["gene"].append(o)
            url = libgly.get_xref_url(map_dict, "protein_xref_uniprotkb", uniprotkb_ac,is_cited)
            o = { "name":species_obj[species]["long_name"],"taxid":species_obj[species]["tax_id"]}
            for nm in ["common_name", "glygen_name"]:
                name_val = species_obj[str(species_obj[species]["tax_id"])][nm]
                if name_val != "":
                    o[nm] = name_val
            o["reference_species"] = "%s [%s]" % (o["name"], o["taxid"])
            protein_obj_dict[canon]["species"].append(o)
    

    file_list = glob.glob("reviewed/*_protein_*names.csv")
    file_list += glob.glob("reviewed/*_protein_*names_refseq.csv")
    file_list += glob.glob("reviewed/*_protein_*names_uniprotkb.csv")

    sheet_info = load_name_sheet_info()
    for in_file in file_list:
        sheet_name = "protein_" + in_file.split("_protein_")[1].replace(".csv", "")
        prop_name = "gene_names" if sheet_name.find("genenames") != -1 else "protein_names"
        field_list = sheet_info[sheet_name]["fieldlist"]
        resource = sheet_info[sheet_name]["resource"]
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            for f in field_list:
                name = row[f_list.index(f)]
                if name.strip() == "":
                    continue
                name_type = "synonym"
                if f in ["recommended_name_full", "gene_symbol_recommended"]:
                    name_type = "recommended"
                #if f == "recommended_name_full":
                #    recname_dict["protein"][canon] = name
                #if f == "gene_symbol_recommended":
                #    recname_dict["gene"][canon] = name
                #if sheet_name == "protein_submittednames":
                #    submittedname_dict["protein"][canon] = name
                xref_key = sheet_info[sheet_name]["xref_key"]
                xref_id = canon.split("-")[0]
                if xref_key == "protein_xref_refseq":
                    xref_id = row[f_list.index("refseq_ac")]
                xref_badge = libgly.get_xref_badge(map_dict, xref_key) 
                xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                name_obj = { "name":name, "resource":xref_badge,  "type":name_type,
                        "url":xref_url, "id":xref_id}

                if canon not in protein_obj_dict:
                    protein_obj_dict[canon] = {"gene":[], "species":[]}
                if prop_name not in protein_obj_dict[canon]:
                    protein_obj_dict[canon][prop_name] = []
                protein_obj_dict[canon][prop_name].append(name_obj)

    return


def load_name_sheet_info():

    return {
        "protein_recnames":{
            "resource":"uniprotkb",
            "xref_key":"protein_xref_uniprotkb_proteinname",
            "fieldlist":["recommended_name_full","recommended_name_short","ec_name"]
        },
        "protein_altnames":{
            "resource":"uniprotkb",
            "xref_key":"protein_xref_uniprotkb_proteinname",
            "fieldlist":["alternative_name_full","alternative_name_short","ec_name"]
        },
        "protein_submittednames":{
            "resource":"uniprotkb",
            "xref_key":"protein_xref_uniprotkb_proteinname",
            "fieldlist":["submitted_name_full","submitted_name_short","ec_name"]
        },
        "protein_proteinnames_refseq":{
            "resource":"refseq",
            "xref_key":"protein_xref_refseq",
            "fieldlist":["refseq_protein_name"]
        },
        "protein_genenames_uniprotkb":{
            "resource":"uniprotkb",
            "xref_key":"protein_xref_uniprotkb",
            "fieldlist":["gene_symbol_recommended","gene_symbol_alternative","orf_name"]
        },
        "protein_genenames_refseq":{
            "resource":"refseq",
            "xref_key":"protein_xref_refseq",
            "fieldlist":["refseq_gene_name"]
        }
    }



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
    sec_name_list = []
    if options.sec != None:
        sec_name_list = options.sec.split(",")

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global main_dict
    global is_cited



    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"

    DEBUG = False




    log_file = "logs/make-glycandb.log"
    msg = "make-glycandb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")


    species_obj, species_list = {}, []
    load_species_info(species_obj, species_list)

    motif_name_dict = load_motif_names()



    motif2gtc = load_motif2gtc()
    motif_glytoucan_ac_list = list(motif2gtc.keys())


    is_cited = libgly.get_is_cited()

    #Clear existing jumbo objects of this record type
    cmd = "rm -f jsondb/jumbodb/glycandb/*"
    x = subprocess.getoutput(cmd)
    cmd = "rm -f jsondb/batchdb/glycan.*"
    x = subprocess.getoutput(cmd)




    sec_info = json.loads(open("generated/misc/glycan_sectioninfo.json", "r").read())
    main_dict = {}
    for sec in sec_info:
        main_dict[sec] = {"seen":{}}

    map_dict = {}
    csvutil.load_dictionaries(map_dict, "generated/misc/")


    residue_heirarchy = json.loads(open("generated/misc/residue_heirarchy.json", "r").read())
    residue2class = {}
    for c in residue_heirarchy:
        for m in residue_heirarchy[c]:
            residue2class[m] = c

    fully_determined_list = []
    data_frame = {}
    in_file = "reviewed/glycan_fully_determined.csv" 
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        fully_determined_list.append(main_id)
        
    
    # set order to 1000 to move datasets to the end
    # set order to 0 to move datasets to the start
    order_dict = {"glycan_masterlist":0}
    expected_dslist =  csvutil.get_expected_dslist("glycan", sec_name_list, order_dict)
    
    #print (json.dumps(expected_dslist, indent=4))
    #print (len(expected_dslist))
    #exit()





   
    biomarker_dict = csvutil.get_biomarker_dict("glycan")

    
    #in_file = "%s%s" % (data_dir,"glycan_subsumption.csv")
    #heirarchy_dict = load_subsumption_heirarchy(in_file)

    protein_obj_dict = {}
    if DEBUG == False:
        load_protein_objects(protein_obj_dict)



    record_stat = {}
   

    missing_list = []
    for file_name in expected_dslist:
        in_file = "%s%s.csv" % (data_dir,file_name)
        if os.path.isfile(in_file) == False:
            missing_list.append(file_name)
    if missing_list != []:
        print ("Dataset files missing:")
        print (json.dumps(missing_list, indent=4))
        sys.exit()

    #expected_dslist = list(set(expected_dslist) - set(missing_list))








    file_idx = 1
    file_count = len(expected_dslist)
    main_obj_dict = {}
    for file_name in expected_dslist:
        in_file = "%s%s.csv" % (data_dir,file_name)
        msg = "make-glycandb: %s [%s/%s]" % (in_file, file_idx, file_count)
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
                missing_score = tmp_row[tmp_fl.index("missing_score")]
                number_monosaccharides_one = tmp_row[tmp_fl.index("monosaccharides")]
                glycan_type = tmp_row[tmp_fl.index("glytoucan_type")]
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
                number_monosaccharides_two = number_monosaccharides_one.replace("+", "")
                if number_monosaccharides_two.isdigit():
                    prop_name = "number_monosaccharides"
                    main_dict[prop_name][main_id] = int(number_monosaccharides_two)
                    combo_list = ["glytoucan_ac", "number_monosaccharides"]
                    update_record_stat(record_stat, file_name, prop_name, 1,combo_list)
                    if number_monosaccharides_one.find("+") != -1:
                        prop_name = "number_monosaccharides_suffix"
                        tmp_val = number_monosaccharides_one.replace(number_monosaccharides_two,"")
                        main_dict[prop_name][main_id] = tmp_val
                        combo_list = ["glytoucan_ac", "number_monosaccharides_suffix"]
                        update_record_stat(record_stat, file_name, prop_name, 1,combo_list)

                if glycan_type != "":
                    prop_name = "glycan_type"
                    main_dict[prop_name][main_id] = glycan_type
                    combo_list = ["glytoucan_ac", "glycan_type"]
                    update_record_stat(record_stat, file_name, prop_name, 1,combo_list)
                if missing_score != "":
                    prop_name = "missing_score"
                    main_dict[prop_name][main_id] = int(missing_score)
                    combo_list = ["glytoucan_ac", "missing_score"]
                    update_record_stat(record_stat, file_name, prop_name, 1,combo_list)
                    if main_dict[prop_name][main_id] == 0:
                        kw = "fully_determined"
                        if main_id not in main_dict["keywords"]:
                            main_dict["keywords"][main_id] = []
                        if kw not in main_dict["keywords"][main_id]:
                            main_dict["fully_determined"][main_id] = "yes"
                            main_dict["keywords"][main_id].append(kw)

        #--> fully_determined
        sheet_name = "fully_determined"
        if file_name.find(sheet_name) != -1:
            species = file_name.split("_")[0]
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
            tmp_fl = data_frame["fields"]
            for main_id in data_frame["data"]:
                main_dict["fully_determined"][main_id] = "yes"
                kw = "fully_determined"
                if main_id not in main_dict["keywords"]:
                    main_dict["keywords"][main_id] = []
                if kw not in main_dict["keywords"][main_id]:
                    main_dict["keywords"][main_id].append(kw)


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
                "glycam":"sequence_glycam_iupac",
                "byonic":"sequence_byonic",
                "gwb":"sequence_gwb"
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
                            url = libgly.get_xref_url(map_dict, "glycan_xref_inchi_key", inchi_key,is_cited)
                            o = {"key":inchi_key, "url":url}
                            main_dict["inchi_key"][main_id] = o
                            combo_list = ["glytoucan_ac", field]
                            update_record_stat(record_stat, file_name, "inchi_key", 1, combo_list)
        
        
        #--> species
        for sheet_name in ["glycan_species"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "species"
                prop_dict = {"taxid":"tax_id", "name":"tax_name", "annotation_category":"annotation_category"}
                xref_info = ["xref_key", "xref_id"]
                combo_flist_one = ["tax_id", "annotation_category"]
                combo_flist_two = ["xref_key", "xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj,in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two, "glytoucan_ac")
                combo_list = ["glytoucan_ac"] + combo_flist_one
                update_record_stat(record_stat, file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id, tax_id,ann_cat = combo_id.split("|")
                    if load_obj[combo_id]["annotation_category"] not in ["Direct","Subsumption"]:
                        load_obj[combo_id]["deleteflag"] = True
                    
                    #Consider only reference species
                    #if species_obj[tax_id]["is_reference"] == "no":
                    #    load_obj[combo_id]["deleteflag"] = True
                    #load_obj[combo_id].pop("annotation_category")
                    if load_obj[combo_id]["taxid"]  == "":
                        load_obj[combo_id]["taxid"] = -1
                    else:
                        tax_id = load_obj[combo_id]["taxid"]
                        long_name = species_obj[str(tax_id)]["long_name"]
                        load_obj[combo_id]["name"] = long_name
                        load_obj[combo_id]["taxid"] = int(tax_id)
                        for nm in ["common_name", "glygen_name"]:
                            name_val = species_obj[str(tax_id)][nm]
                            load_obj[combo_id][nm] = name_val
                        reference_species = "%s [%s]" % (long_name, tax_id)
                        load_obj[combo_id]["reference_species"] = reference_species
           



                    


        #--> interactions
        for sheet_name in ["protein_matrixdb"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "interactions"
                prop_dict = {"interactor_id":"uniprotkb_canonical_ac"}
                xref_info = ["xref_key", "xref_id"]
                combo_flist_one = ["uniprotkb_canonical_ac"]
                combo_flist_two = ["xref_key", "xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj,in_file, prop_dict,xref_info,
                        combo_flist_one, combo_flist_two, "saccharide")
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    if combo_id.split("|")[1] == "":
                        load_obj[combo_id]["deleteflag"] = True

                    canon = combo_id.split("|")[1]
                    rec_name = ""
                    if canon in protein_obj_dict:
                        protein_obj = protein_obj_dict[canon]
                        rec_name = get_protein_info(protein_obj, "recommendedname")
                    load_obj[combo_id]["interactor_type"] = "protein"
                    load_obj[combo_id]["interactor_name"] = rec_name

        #--> publication
        for sheet_name in ["glycan_citations_", "_proteoform_citations_"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "publication"
                prop_dict = {"title":"title", 
                        "journal":"journal_name","date":"publication_date","authors":"authors"}
                xref_info = ["src_xref_key", "src_xref_id"]
                #combo_flist_one = ["xref_key", "xref_id"]
                combo_flist_one = ["xref_id"]
                combo_flist_two = ["src_xref_key", "src_xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj,in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two, "glytoucan_ac") 
                combo_list = ["glytoucan_ac"] + combo_flist_one
                update_record_stat(record_stat, file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    #main_id, xref_key, xref_id = combo_id.split("|")
                    main_id, xref_id = combo_id.split("|")
                    xref_key = "glycan_xref_pubmed" if xref_id.isdigit() else "glycan_xref_doi" 
                    #print ("Robel-A", combo_id, file_name)
                    xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    load_obj[combo_id]["date"] = load_obj[combo_id]["date"].split(" ")[0]
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    load_obj[combo_id]["reference"] = [
                        {
                            "type":xref_badge,
                            "id":xref_id,
                            "url":xref_url
                        }
                    ] 


        #--> motifs
        for sheet_name in ["glycan_motif"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "motifs"
                load_obj = main_dict[prop_name]
                prop_dict = {"id":"motif_ac", "name":"motif_name", "synonym":"alternative_name",
                        "keywords":"keyword"}
                combo_flist_one = ["motif_ac"]
                xref_info = []
                combo_flist_two = []
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                        combo_flist_one, combo_flist_two, "glytoucan_ac")
                combo_list = ["glytoucan_ac"] + combo_flist_one
                update_record_stat(record_stat, file_name, prop_name, n, combo_list)



        #--> glycoprotein
        for sheet_name in ["proteoform_glycosylation_sites"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "glycoprotein"
                load_obj = main_dict[prop_name]
                prop_dict = {"uniprot_canonical_ac":"uniprotkb_canonical_ac",
                    "start_pos":"glycosylation_site_uniprotkb",
                    "end_pos":"glycosylation_site_uniprotkb",
                    "residue":"amino_acid"
                }
                combo_flist_one = ["uniprotkb_canonical_ac", "glycosylation_site_uniprotkb"]
                #xref_info = ["xref_key", "xref_id"]
                xref_info = ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
                combo_flist_two = ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
                if file_name.find("proteoform_glycosylation_sites_gptwiki") != -1:
                    xref_info = ["glycan_xref_key", "glycan_xref_id"]
                    combo_flist_two = ["glycan_xref_key", "glycan_xref_id"]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two, "saccharide")
                combo_list = ["glytoucan_ac"] + combo_flist_one
                update_record_stat(record_stat, file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id, canon, aa_pos = combo_id.split("|")[:3]
                    if aa_pos.strip() == "":
                        if "start_pos" in load_obj[combo_id]:
                            load_obj[combo_id].pop("start_pos")
                            load_obj[combo_id].pop("end_pos")
                            load_obj[combo_id].pop("residue") 
                    protein_obj = protein_obj_dict[canon] if canon in protein_obj_dict else {}
                    if protein_obj == {}:
                        load_obj[combo_id]["deleteflag"] = True
                        continue
                    gene_name = get_protein_info(protein_obj, "gene_name")
                    rec_name = get_protein_info(protein_obj, "recommendedname")
                    tax_id = get_protein_info(protein_obj, "tax_id")
                    tax_name = get_protein_info(protein_obj, "tax_name")
                    tax_common_name = get_protein_info(protein_obj, "tax_common_name")
                    load_obj[combo_id]["protein_name"] = rec_name
                    load_obj[combo_id]["gene_name"] = gene_name
                    load_obj[combo_id]["tax_id"] = tax_id
                    load_obj[combo_id]["tax_name"] = tax_name
                    load_obj[combo_id]["tax_common_name"] = tax_common_name
                    if "start_pos" in load_obj[combo_id]:
                        s_pos_list = str(load_obj[combo_id]["start_pos"]).split("|")
                        e_pos_list = str(load_obj[combo_id]["end_pos"]).split("|")
                        load_obj[combo_id]["start_pos"] = int(s_pos_list[0])
                        load_obj[combo_id]["end_pos"] = int(e_pos_list[0])
                        if len(s_pos_list) > 1:
                            load_obj[combo_id]["alternate_start_pos_list"] = []
                            for p in s_pos_list[1:]:
                                load_obj[combo_id]["alternate_start_pos_list"].append(int(p))
                        if len(e_pos_list) > 1:
                            load_obj[combo_id]["alternate_end_pos_list"] = []
                            for p in e_pos_list[1:]:
                                load_obj[combo_id]["alternate_end_pos_list"].append(int(p))
                    if "residue" in load_obj[combo_id]:
                        r_list = load_obj[combo_id]["residue"].split("|")
                        load_obj[combo_id]["residue"] = r_list[0]
                        if len(r_list) > 1:
                            load_obj[combo_id]["alternate_residue_list"] = s_pos_list[1:]
                                                                                
        

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
                        #if domain.lower() == "byonic":
                        #    main_dict["byonic"][main_id] = name



        #--> subsumption
        for sheet_name in ["glycan_subsumption"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "subsumption"
                rs_dict = {
                    "ancestor":1,"basecomposition":1, "composition":1,"descendant":1,
                    "fullydetermined":1,"subsumedby":1,"subsumes":1,"topology":1,"leaf":1
                }
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        related_accession = tmp_row[tmp_fl.index("related_accession")]
                        relationship = tmp_row[tmp_fl.index("relationship")].lower()
                        if relationship in rs_dict:
                            combo_id = "%s|%s|%s" % (main_id, related_accession,relationship)
                            o = {"related_accession":related_accession,"relationship":relationship}
                            load_obj[combo_id] = o

                        #glytoucan_type = tmp_row[tmp_fl.index("glytoucan_type")].lower()
                        #gt_list = ["topology", "composition", "basecomposition"]
                        #rl_list = ["ancestor", "descendant"]
                        #cond_one = main_id != related_accession
                        #cond_two = glytoucan_type == relationship
                        #cond_three = glytoucan_type in gt_list
                        #if cond_one and cond_two and cond_three:
                        #    combo_id = "%s|%s|%s" % (main_id, related_accession,relationship)
                        #    rl = ""
                        #    if main_id in heirarchy_dict:
                        #        if related_accession in heirarchy_dict[main_id]:
                        #            rl = heirarchy_dict[main_id][related_accession]
                        #    load_obj[combo_id] = {
                        #        "id":related_accession, 
                        #        "type":glytoucan_type, 
                        #        "relationship":rl
                        #    }


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


        #--> enzyme
        for sheet_name in ["glycan_enzyme"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "enzyme"
                load_obj = main_dict[prop_name]
                prop_dict = {"uniprot_canonical_ac":"uniprotkb_canonical_ac", 
                        "xref_key":"xref_key", "xref_id":"xref_id"}
                combo_flist_one = ["uniprotkb_canonical_ac"]
                xref_info = []
                combo_flist_two = []
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two, "glytoucan_ac")
                combo_list = ["glytoucan_ac"] + combo_flist_one
                update_record_stat(record_stat, file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id, canon = combo_id.split("|")
                    protein_obj = protein_obj_dict[canon] if canon in protein_obj_dict else {}
                    if protein_obj == {}:
                        continue
                    gene_name = get_protein_info(protein_obj, "gene_name")
                    rec_name = get_protein_info(protein_obj, "recommendedname")
                    gene_url = get_protein_info(protein_obj, "gene_url")
                    enzyme_tax_id = get_protein_info(protein_obj, "tax_id")
                    enzyme_tax_name = get_protein_info(protein_obj, "tax_name")
                    enzyme_tax_common_name = get_protein_info(protein_obj, "tax_common_name")
                    load_obj[combo_id]["protein_name"] = rec_name
                    load_obj[combo_id]["gene"] = gene_name
                    load_obj[combo_id]["gene_link"] = gene_url
                    load_obj[combo_id]["tax_id"] = enzyme_tax_id
                    load_obj[combo_id]["tax_name"] = enzyme_tax_name
                    load_obj[combo_id]["tax_common_name"] = enzyme_tax_common_name
                    xref_key = load_obj[combo_id]["xref_key"]
                    xref_id = load_obj[combo_id]["xref_id"]
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    ev_obj = {"id":xref_id, "database":xref_badge, "url":xref_url}
                    load_obj[combo_id]["evidence"] = [ev_obj]


        #--> tool_support
        for sheet_name in ["glycan_toolsupport"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "tool_support"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        tool = tmp_row[tmp_fl.index("tool")]
                        support = tmp_row[tmp_fl.index("support")]
                        if main_id not in load_obj:
                            load_obj[main_id] = {}
                        load_obj[main_id][tool] = support


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

        #--> glycotree_pathways
        for sheet_name in ["glycan_pathway_glycotree"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "glycotree_pathways"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        val_dict = {}
                        for f in tmp_fl:
                            val_dict[f] = tmp_row[tmp_fl.index(f)]
                        o = {
                            "source":val_dict["source"],
                            "target":val_dict["target"],
                            "residue_affected": {
                                "id": val_dict["residue_id"],
                                "full_name": val_dict["residue_name"]
                            },
                            "enzymes":[]
                        }
                        combo_id = "%s|%s|%s" % (main_id, val_dict["source"], val_dict["target"])
                        load_obj[combo_id] = o
                        ac, tax_name = val_dict["enzyme_uniprotkb_ac"], val_dict["enzyme_tax_name"]
                        if ac != "":
                            oo = {"uniprotkb_ac":ac,"tax_name":tax_name}
                            load_obj[combo_id]["enzymes"].append(oo)



        #--> crossref
        for sheet_name in ["glycan_xref_"]:
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
                    xref_key, xref_id = combo_id.split("|")[-2], combo_id.split("|")[-1]
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    load_obj[combo_id]["url"] = xref_url
                    load_obj[combo_id]["database"] = xref_badge
                    cats_dict = map_dict["xrefkey2category"]
                    xref_categories = cats_dict[xref_key] if xref_key in cats_dict else []
                    load_obj[combo_id]["categories"] = xref_categories

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
                            
        #--> expression
        for sheet_name in ["proteoform_glycosylation_sites_"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "expression"
                load_obj = main_dict[prop_name]
                FR = open(in_file, "r")
                idx = 0
                tmp_fl = []
                for line in FR:
                    idx += 1
                    tmp_row = line.strip().split("\",\"")
                    tmp_row[0], tmp_row[-1] = tmp_row[0].replace("\"", ""), tmp_row[-1].replace("\"", "")
                    if idx == 1:
                        tmp_fl = tmp_row
                        continue
                    else:
                        main_id = tmp_row[tmp_fl.index("saccharide")]
                        if main_id == "":
                            continue
                        canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                        tmp_aa_pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                        tmp_aa_three = tmp_row[tmp_fl.index("amino_acid")]
                        aa_pos_list = tmp_aa_pos.split("|")
                        aa_three_list = tmp_aa_three.split("|")
                        aa_pos = aa_pos_list[0]
                        aa_three = aa_three_list[0]
                        #if aa_pos == "":
                        #    continue
                        abundance, tissue_name, tissue_id, tissue_namespace = "", "", "", ""
                        if "source_tissue_id" in tmp_fl:
                            tissue_name = tmp_row[tmp_fl.index("source_tissue_name")]
                            t = tmp_row[tmp_fl.index("source_tissue_id")]
                            if t.find(":") != -1:
                                tissue_namespace, tissue_id = t.split(":")[0], t.split(":")[1]
                            elif t.find("_") != -1:
                                tissue_namespace, tissue_id = t.split("_")[0], t.split("_")[1]
       
                        if "abundance" in tmp_fl:
                            abundance = tmp_row[tmp_fl.index("abundance")]
                        
                        cl_name, cl_id, cl_namespace = "", "", ""
                        if "source_cell_line_cellosaurus_id" in tmp_fl:
                            cl_name = tmp_row[tmp_fl.index("source_cell_line_cellosaurus_name")]
                            c = tmp_row[tmp_fl.index("source_cell_line_cellosaurus_id")]
                            if c.find(":") != -1:
                                cl_namespace, cl_id = c.split(":")[0], c.split(":")[1]
                            elif c.find("_") != -1:
                                cl_namespace, cl_id = c.split("_")[0], c.split("_")[1]

                        src_xref_key = tmp_row[tmp_fl.index("src_xref_key")]
                        src_xref_id = tmp_row[tmp_fl.index("src_xref_id")]
                        if src_xref_key == "protein_xref_glyconnect":
                            src_xref_key = "glycan_xref_glyconnect"
                            src_xref_id = tmp_row[tmp_fl.index("structure_id")]
                        src_xref_badge = libgly.get_xref_badge(map_dict, src_xref_key)
                        src_xref_url =  libgly.get_xref_url(map_dict, src_xref_key, src_xref_id,is_cited)
                        xref_key = tmp_row[tmp_fl.index("xref_key")]
                        xref_id = tmp_row[tmp_fl.index("xref_id")]
                        xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                       
                        cl_obj, tissue_obj = {}, {}
                        if cl_id not in ["", "0"]:
                            cl_key = "cell_xref_" + cl_namespace.lower()
                            if cl_key in map_dict["xrefkey2url"]:
                                cl_url =  libgly.get_xref_url(map_dict, cl_key, cl_id,is_cited)
                                cl_obj = {"name": cl_name,"namespace": cl_namespace,"id":cl_id,"url":cl_url}
                        if tissue_id != "":
                            t_ns, t_nm = tissue_namespace, tissue_name
                            t_xref_key = "tissue_xref_" + tissue_namespace.lower()
                            if t_xref_key in map_dict["xrefkey2url"]:
                                t_xref_url =  libgly.get_xref_url(map_dict, t_xref_key, tissue_id,is_cited)
                                tissue_obj = {"name":t_nm,"namespace":t_ns,"id":tissue_id,"url":t_xref_url}
                       

                         
                        if cl_obj == {} and tissue_obj == {}:
                            continue
                        # this means cell line was derived from tissue, so ignore tissue
                        if cl_obj != {} and tissue_obj != {}:
                            tissue_obj = {}

                        sample_src_id = cl_id if tissue_obj == {} else tissue_id                        
 
                        #Making sure only on of these objects can have values
                        category = "" 
                        if tissue_obj != {}:
                            cl_obj = {}
                            category = "tissue"
                        if cl_obj != {}:
                            tissue_obj = {}
                            category = "cell_line"
                        
                        combo_id = "%s|%s|%s|%s|%s" % (main_id,canon,aa_pos,sample_src_id,abundance)
                        if combo_id not in load_obj["seen"]:
                            load_obj[combo_id] = {
                                "category":category,
                                "tissue":tissue_obj,
                                "cell_line":cl_obj,
                                "abundance":abundance,
                                "evidence":[]
                            }
                            if aa_pos != "":
                                load_obj[combo_id]["start_pos"], load_obj[combo_id]["end_pos"] = int(aa_pos), int(aa_pos)
                                load_obj[combo_id]["residue"] = aa_three
                            if canon != "":
                                load_obj[combo_id]["uniprot_canonical_ac"] = canon
                            if len(aa_pos_list) > 1:
                                load_obj[combo_id]["alternate_start_pos_list"] = []
                                load_obj[combo_id]["alternate_end_pos_list"] = []
                                for tmp_aa_pos in aa_pos_list[1:]:
                                    p = int(tmp_aa_pos)
                                    load_obj[combo_id]["alternate_start_pos_list"].append(p)
                                    load_obj[combo_id]["alternate_end_pos_list"].append(p)
                            if len(aa_three_list) > 1:
                                load_obj[combo_id]["alternate_residue"] = aa_three_list[1:]
                            load_obj["seen"][combo_id] = True

                        evdn_row_list = [ [src_xref_key,src_xref_id], [xref_key,xref_id]]
                        for evdn_row in evdn_row_list:
                            x_key, x_id = evdn_row[0], evdn_row[1]
                            x_badge = libgly.get_xref_badge(map_dict, x_key)
                            x_url =  libgly.get_xref_url(map_dict, x_key, x_id,is_cited)
                            if x_id == "" or x_url == "":
                                continue
                            evdn_combo = "%s|%s|%s|%s|%s|%s|%s"%(main_id,canon,aa_pos,sample_src_id,abundance,x_key,x_id)
                            if evdn_combo not in load_obj["seen"]:
                                ev_obj = {"id":x_id, "database":x_badge, "url":x_url}
                                load_obj[combo_id]["evidence"].append(ev_obj)
                                load_obj["seen"][evdn_combo] = True
                FR.close()

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
                                url =  libgly.get_xref_url(map_dict, "glycan_xref_monosaccharide_residue_name",cid,is_cited)
                                o = {"name":name,"residue":residue,"count":n,"cid":cid,"url":url}

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
                            if residue != "Sia" and residue not in residue_heirarchy:
                                continue
                            r_class = residue
                            n = int(tmp_row[tmp_fl.index(residue)])
                            class_count[r_class] = n
                       
                        n_sia = class_count["Sia"] if "Sia" in class_count else 0
                        n_neuac = class_count["NeuAc"] if "NeuAc" in class_count else 0
                        n_neugc = class_count["NeuGc"] if "NeuGc" in class_count else 0
                        if n_sia - n_neuac - n_neugc > 0:
                            class_count["Xxx"] += n_sia - n_neuac - n_neugc
                        class_count.pop("Sia")
 
                        for r_class in class_count:
                            n = class_count[r_class]
                            r_class_name = map_dict["residue2name"][r_class][0]
                            cid = map_dict["residue2cid"][r_class][0]
                            r_class = "other" if r_class.lower() == "xxx" else r_class.lower()
                            o = {"name":r_class_name, "residue":r_class, "count":n}
                            if cid != "" and cid.strip() != "0":
                                url = libgly.get_xref_url(map_dict, "glycan_xref_monosaccharide_residue_name", cid,is_cited)
                                o = {"name":r_class_name, "residue":r_class, "count":n, "cid":cid, "url":url}
                            combo_id = "%s|%s|%s" % (main_id, r_class, cid)
                            load_obj[combo_id] = o
                combo_list = ["glytoucan_ac", "residue"]
                update_record_stat(record_stat, file_name, prop_name, len(load_obj.keys()), combo_list)

        
        #--> classification
        for sheet_name in ["glycan_classification"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "classification"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    type_combo_list = []
                    for tmp_row in data_frame["data"][main_id]:
                        g_type = tmp_row[tmp_fl.index("glycan_type")].strip()
                        g_subtype = tmp_row[tmp_fl.index("glycan_subtype")].strip()
                        g_type = "Other" if g_type == "" else g_type
                        g_subtype = "Other" if g_subtype == "" else g_subtype
                        type_combo_list.append("%s|%s" % (g_type, g_subtype))
                    type_combo_list = sorted(set(type_combo_list))
                    if len(type_combo_list) == 1:
                        glycan_type, glycan_subtype = type_combo_list[0].split("|")
                        t_tag = map_dict["subtypetags"][glycan_type][0] if glycan_type in map_dict["subtypetags"] else "xxx"
                        #type_url = libgly.get_xref_url(map_dict, "glycan_xref_glycan_type", t_tag,is_cited)
                        type_url = "" 
                        o = {"type":{"name":glycan_type, "url":type_url}}
                        if t_tag == "xxx":
                            o["type"].pop("url")
                        if glycan_subtype != "no subtype":
                            s_tag = map_dict["subtypetags"][glycan_subtype][0] if glycan_subtype in map_dict["subtypetags"] else "xxx"
                            #subtype_url = libgly.get_xref_url(map_dict, "glycan_xref_glycan_type", s_tag,is_cited)
                            subtype_url = ""
                            o["subtype"] = {"name":glycan_subtype, "url":subtype_url}
                            if s_tag == "xxx":
                                o["subtype"].pop("url")
                        combo_id = "%s|%s|%s" % (main_id, glycan_type, glycan_subtype)
                        load_obj[combo_id] = o
                    else:
                        for type_combo in type_combo_list:
                            glycan_type, glycan_subtype = type_combo.split("|")
                            if glycan_subtype == "no subtype":
                                continue
                            t_tag = map_dict["subtypetags"][glycan_type][0] if glycan_type in map_dict["subtypetags"] else "xxx"
                            #type_url = libgly.get_xref_url(map_dict, "glycan_xref_glycan_type", t_tag,is_cited)
                            type_url = ""
                            s_tag = map_dict["subtypetags"][glycan_subtype][0] if glycan_subtype in map_dict["subtypetags"] else "xxx"
                            #subtype_url = libgly.get_xref_url(map_dict,  "glycan_xref_glycan_type", s_tag,is_cited)
                            subtype_url = ""
                            o = {
                                "type":{"name":glycan_type, "url":type_url}
                                ,"subtype":{"name":glycan_subtype, "url":subtype_url}
                            }
                            if t_tag == "xxx":
                                o["type"].pop("url")
                            if s_tag == "xxx":
                                o["subtype"].pop("url")
                            combo_id = "%s|%s|%s" % (main_id, glycan_type, glycan_subtype)
                            load_obj[combo_id] = o
                combo_list = ["glytoucan_ac", "glycan_type", "glycan_subtype"]
                update_record_stat(record_stat, file_name, prop_name, len(load_obj.keys()), combo_list)


    #prop_name = "publication"
    #for combo_id in main_dict[prop_name]:
    #    print ("Robel-B", combo_id)



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
        #if main_id in motif_glytoucan_ac_list:
        #    continue

        obj = tmp_dict[main_id]
        if obj["classification"] == []:
            o = {"type":{"name":"Other", "url":""}, "subtype":{"name":"Other", "url":""}}
            obj["classification"] = [o]

        if obj["species"] != []:
            sp_dict = {}
            for o in obj["species"]:
                if o["annotation_category"] == "Subsumption":
                    for jj in range(0, len(o["evidence"])):
                        o["evidence"][jj]["database"] = o["evidence"][jj]["database"].replace("GNOme", "Subsumption") 

                tax_id = o["taxid"]
                if tax_id not in sp_dict:
                    sp_dict[tax_id] = o
                else:
                    sp_dict[tax_id]["evidence"] += o["evidence"]
                
            obj["species"] = []
            for tax_id in sp_dict:
                seen_evdn_db = {}
                for jj in range(0, len(sp_dict[tax_id]["evidence"])):
                    seen_evdn_db[sp_dict[tax_id]["evidence"][jj]["database"]] = True
                cat_str = ";".join(list(seen_evdn_db.keys()))
                if cat_str == "Subsumption":
                    sp_dict[tax_id]["annotation_category"] = "subsumption"
                elif cat_str.find("Subsumption") == -1:
                    sp_dict[tax_id]["annotation_category"] = "direct"
                else:
                    sp_dict[tax_id]["annotation_category"] = "both"
                obj["species"].append(sp_dict[tax_id])

        for k in ["mass", "mass_pme", "number_monosaccharides"]:
            if obj[k] == 0.0 or obj[k] == 0:
                obj.pop(k)
        
        section_stats.get_sec_stats(obj, "glycan")

        # if main_id is motif, at motif name to names
        if main_id in motif_name_dict:
            motif_name = motif_name_dict[main_id]["name"]
            if motif_name.strip() != "":
                obj["names"] += [{"name":motif_name, "domain":"motifname"}]
            for motif_syn in motif_name_dict[main_id]["synonym_list"]:
                if motif_syn.strip() != "":
                    obj["names"] += [{"name":motif_syn, "domain":"motifsynonym"}]

        out_file = "jsondb/glycandb/%s.json" % (main_id)
        out_str = json.dumps(obj, indent=4)
        if len(out_str) > 16000000:
            out_file = "jsondb/jumbodb/glycandb/%s.json" % (main_id)
        
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        record_count += 1 
    msg = "make-glycandb: final filtered in: %s glycan objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "a")
   
    

   




if __name__ == '__main__':
        main()



