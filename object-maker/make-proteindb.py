#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
import subprocess
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq


import libgly
import csvutil
import batchutil
import section_stats






def get_seq_features(obj):

    tmp_dict = {}

    sec = "glycosylation"
    tmp_dict["n_linked_sites"] = []
    tmp_dict["o_linked_sites"] = []
    if sec in obj:
        for o in obj[sec]:
            if "start_pos" not in o:
                continue
            if o["type"].lower() == "n-linked":
                if o["start_pos"] not in  tmp_dict["n_linked_sites"]:
                    tmp_dict["n_linked_sites"].append(o["start_pos"])
            elif o["type"].lower() == "o-linked":
                if o["start_pos"] not in  tmp_dict["o_linked_sites"]:
                    tmp_dict["o_linked_sites"].append(o["start_pos"])
    
    for k in ["n_linked_sites", "o_linked_sites"]:
        tmp_dict[k] = sorted(tmp_dict[k] )


    sec = "glycation"
    for sec in ["snv", "glycation", "phosphorylation"]:
        k = "%s_sites" % (sec)
        tmp_dict[k] = []    
        if sec in obj:
            for o in obj[sec]:
                if "start_pos" not in o:
                    continue
                if o["start_pos"] not in  tmp_dict[k]:
                    tmp_dict[k].append(o["start_pos"])
        tmp_dict[k] = sorted(tmp_dict[k] )

    sec = "site_annotation"
    k = "sequon_annotation_sites"
    tmp_dict[k] = []
    if sec in obj:
        for o in obj[sec]:
            if o["annotation"] in ["n_glycosylation_sequon", "o_glycosylation_sequon"]:
                tmp_dict[k].append({"start_pos":o["start_pos"] , "end_pos":o["end_pos"]})   

    return tmp_dict





##################
def get_sort_key_value_pub(obj):
    return obj["date"]

def get_sort_key_value_mut(obj):
    return obj["ann_score"]

def load_pathways(pathway_dict):

    file_list = glob.glob(data_dir + "/reviewed/*_protein_pathways_reactome.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            p_id = row[f_list.index("pathway_id")]
            reaction_list = []
            if row[f_list.index("reaction_id_list")] != "":
                reaction_list = row[f_list.index("reaction_id_list")].split(",")
            pathway_dict[p_id] = {                   
                "description":row[f_list.index("pathway_name")],
                "reaction_list":reaction_list
            }
    return


def load_enzyme_ann(enzyme_ann, species):
    
    ec_in_rhea = {}
    in_file = data_dir + "/reviewed/protein_reaction2ec_rhea.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ec_number = row[f_list.index("ec_number")]
        ec_in_rhea[ec_number] = True

    ec_in_brenda = {}
    for in_file in glob.glob(data_dir + "/reviewed/*_protein_xref_brenda.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ec_number = row[f_list.index("xref_id")]
            ec_in_brenda[ec_number] = True

    tmp_dict = {}
    file_list = glob.glob(data_dir + "/reviewed/%s_protein_enzyme_annotation_uniprotkb.csv" % (species))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            ec_number = row[f_list.index("enzyme_ec")]
            ec_name = row[f_list.index("enzyme_activity")]
            if ec_number == "" or ec_name == "":
                continue
            evdn_obj_list = []
            if ec_number in ec_in_rhea:
                xref_key, xref_id = "protein_xref_rhea_enzyme", ec_number
                xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                evdn_obj_list.append({"id":xref_id, "database":xref_badge, "url":xref_url})
            if ec_number in ec_in_brenda:
                xref_key, xref_id = "protein_xref_brenda_enzyme", ec_number
                xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                evdn_obj_list.append({"id":xref_id, "database":xref_badge, "url":xref_url})

            if canon not in tmp_dict:
                tmp_dict[canon] = {}
            if ec_number not in tmp_dict[canon]:
                tmp_dict[canon][ec_number] = {"ec_number":ec_number, "ec_name":ec_name, "evidence":evdn_obj_list}
                
    for canon in tmp_dict:
        enzyme_ann[canon] = []
        for ec_number in tmp_dict[canon]:
            enzyme_ann[canon].append(tmp_dict[canon][ec_number])


    return


def load_enzyme_dict(enzyme_dict, species):

    canon2gene_name = {}
    file_list = glob.glob(data_dir + "/reviewed/%s_protein_masterlist.csv" % (species))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            canon2gene_name[canon] = row[f_list.index("gene_name")]

    file_list = glob.glob(data_dir + "/reviewed/%s_protein_enzyme_annotation_uniprotkb.csv" % (species))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            e_id = canon.split("-")[0]
            xref_key, xref_id = "protein_xref_uniprotkb", e_id
            xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
            xref_badge = libgly.get_xref_badge(map_dict, xref_key)
            if canon not in enzyme_dict:
                enzyme_dict[canon] = {
                    "uniprot_ac":e_id,"gene_name":canon2gene_name[canon],
                    "activity":[]
                }
            o = {"ec_number":row[f_list.index("enzyme_ec")],
                "activity":row[f_list.index("enzyme_activity")],
                "evidence":[{"id":xref_id, "database":xref_badge, "url":xref_url}]
            }
            enzyme_dict[canon]["activity"].append(o)

 
    kw = "enzyme"
    for main_id in enzyme_dict:
        if main_id not in main_dict["keywords"]:
            main_dict["keywords"][main_id] = []
            if kw not in main_dict["keywords"][main_id]:
                main_dict["keywords"][main_id].append(kw)
                combo_list = ["uniprotkb_canonical_ac", "keywords"]



    return


def load_reactions(reaction_dict, participant_dict, ac2rxn, enzyme_dict, species):

    role_dict = {}
    file_list = glob.glob(data_dir + "/reviewed/%s_protein_participants_reactome.csv" % (species))
    file_list += glob.glob(data_dir + "/reviewed/%s_protein_participants_rhea.csv" % (species))

    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            rxn_id = row[f_list.index("reaction_id")]
            participant_id = row[f_list.index("participant_id")]
            role = row[f_list.index("role")]
            if participant_id not in participant_dict:
                participant_dict[participant_id] = {
                    "participant_id":row[f_list.index("participant_id")],
                    "participant_name":row[f_list.index("participant_name")],
                    "crossref":[]
                }
            xref_key = row[f_list.index("xref_type")]
            xref_id = row[f_list.index("xref_id")]
            xref_badge = row[f_list.index("xref_type")]
            xref_url = ""
            o = {"id":xref_id, "database":xref_badge, "url":xref_url}
            participant_dict[participant_id]["crossref"].append(o)
            if rxn_id not in role_dict:
                role_dict[rxn_id] = {"input":[], "output":[], "enzyme":[]}
            if participant_id not in role_dict[rxn_id][role]:
                role_dict[rxn_id][role].append(participant_id)
            if role == "enzyme":
                if participant_id not in ac2rxn:
                    ac2rxn[participant_id] = []
                if rxn_id not in ac2rxn[participant_id]:
                    ac2rxn[participant_id].append(rxn_id)
                if participant_id not in enzyme_dict:
                    role_dict[rxn_id][role].remove(participant_id)




    file_list = glob.glob(data_dir + "/reviewed/%s_protein_reactions_reactome.csv" % (species))
    file_list += glob.glob(data_dir + "/reviewed/%s_protein_reactions_rhea.csv" % (species))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        source = in_file.split("_")[-1].split(".")[0]
        for row in data_frame["data"]:
            rxn_id = row[f_list.index("reaction_id")]
            xref_id = rxn_id
            xref_key = "protein_xref_%s" % (source)
            
            xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
            xref_badge = libgly.get_xref_badge(map_dict, xref_key)
            #"reaction_id","reaction_name","cellular_location","pmid","pathway_id"
            input_list = role_dict[rxn_id]["input"] if rxn_id in role_dict else []
            output_list = role_dict[rxn_id]["output"] if rxn_id in role_dict else []
            enzyme_list = role_dict[rxn_id]["enzyme"] if rxn_id in role_dict else []

            reaction_dict[rxn_id] = {
                "id":rxn_id,
                "name":row[f_list.index("reaction_name")],
                "equation":row[f_list.index("equation")],
                "summary":row[f_list.index("reaction_summary")],
                "cellular_location":row[f_list.index("cellular_location")],
                "input_participant_list":input_list,
                "output_participant_list":output_list,
                "enzyme_list":enzyme_list,
                "evidence":[{"id":rxn_id, "database":xref_badge, "url":xref_url}]
            }


    return

def load_protein_names(recname_dict, submittedname_dict):

    file_list = glob.glob(data_dir + "/reviewed/*_protein_*names.csv")
    file_list += glob.glob(data_dir + "/reviewed/*_protein_*names_refseq.csv")
    file_list += glob.glob(data_dir + "/reviewed/*_protein_*names_uniprotkb.csv")
    sheet_info = load_name_sheet_info()
    for in_file in file_list:
        sheet_name = "protein_" + in_file.split("_protein_")[1].replace(".csv", "")
        prop_name = "gene_names" if sheet_name.find("genenames") != -1 else "protein_names"
        load_obj = main_dict[prop_name]
        field_list = sheet_info[sheet_name]["fieldlist"]
        resource = sheet_info[sheet_name]["resource"]
        resource_badge = sheet_info[sheet_name]["badge"]
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
                if f == "recommended_name_full":
                    recname_dict["protein"][canon] = name
                if f == "gene_symbol_recommended":
                    recname_dict["gene"][canon] = name
                if sheet_name == "protein_submittednames":
                    submittedname_dict["protein"][canon] = name
                xref_key = sheet_info[sheet_name]["xref_key"]
                xref_id = canon.split("-")[0]
                if xref_key == "protein_xref_refseq":
                    xref_id = row[f_list.index("refseq_ac")]
                xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                name_obj = { "name":name, "resource":resource_badge,  "type":name_type,
                        "url":xref_url, "id":xref_id}
                combo_id = "%s|%s|%s|%s" % (canon,name_type,name,resource)
                load_obj[combo_id] = name_obj


    return



def get_sorting_key(obj):
    return obj['sortorder']



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

    #tax_id,short_name,long_name,common_name,glygen_name,nt_file,is_reference,sort_order



    species_obj_list.sort(key=get_sorting_key)
    for o in species_obj_list:
        species_list.append(o["shortname"])


    return




def load_disease_names():

    data_frame = {}
    in_file = data_dir + "/reviewed/protein_disease_names.csv"
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "xref_id")
    tmp_fl = data_frame["fields"]
    for xref_id in data_frame["data"]:
        for tmp_row in data_frame["data"][xref_id]:
            xref_key = tmp_row[tmp_fl.index("xref_key")]
            name_type = tmp_row[tmp_fl.index("name_type")]
            name = tmp_row[tmp_fl.index("name")]
            disease_desc = tmp_row[tmp_fl.index("description")]
            if xref_key not in disease_id2names:
                disease_id2names[xref_key] = {}
                disease_id2desc[xref_key] = {}
            if xref_id not in disease_id2names[xref_key]:
                disease_id2names[xref_key][xref_id] = {"name":"", "synonyms":[]}
            if name_type == "recommended_name":
                disease_id2names[xref_key][xref_id]["name"] = name
                disease_id2desc[xref_key][xref_id] = disease_desc
            else:
                if name not in disease_id2names[xref_key][xref_id]["synonyms"]:
                    disease_id2names[xref_key][xref_id]["synonyms"].append(name)

    return


def load_disease_idmap(doid2xrefid, xrefid2doid):

    data_frame = {}
    in_file = data_dir + "/reviewed/protein_disease_idmap.csv"
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "do_id")
    tmp_fl = data_frame["fields"]
    for do_id in data_frame["data"]:
        for tmp_row in data_frame["data"][do_id]:
            xref_key = tmp_row[tmp_fl.index("xref_key")]
            xref_id = tmp_row[tmp_fl.index("xref_id")]
            if do_id not in doid2xrefid:
                doid2xrefid[do_id] = {}
            if xref_key not in doid2xrefid[do_id]:
                doid2xrefid[do_id][xref_key] = []
            if xref_id not in doid2xrefid[do_id][xref_key]:
                doid2xrefid[do_id][xref_key].append(xref_id)
            if xref_key not in xrefid2doid:
                xrefid2doid[xref_key] = {}
            if xref_id not in xrefid2doid[xref_key]:
                xrefid2doid[xref_key][xref_id] = []
            if do_id not in xrefid2doid[xref_key][xref_id]:
                xrefid2doid[xref_key][xref_id].append(do_id)


    data_frame = {}
    in_file = data_dir + "/reviewed/human_protein_genomics_england_disease.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    tmp_fl = data_frame["fields"]
    xref_key = "protein_xref_genomics_england"
    xrefid2doid[xref_key] = {}
    for row in data_frame["data"]:
        do_id = row[tmp_fl.index("doid")]
        xref_id = row[tmp_fl.index("gene_name")]
        if xref_id not in xrefid2doid[xref_key]:
            xrefid2doid[xref_key][xref_id] = []
        if do_id not in xrefid2doid[xref_key][xref_id]:
            xrefid2doid[xref_key][xref_id].append(do_id)

    return







def load_property_objlist(tmp_obj_dict,in_file, prop_dict,xref_info, combo_flist_one, combo_flist_two):

    record_count = 0
    local_seen_dict = {}
    seen_dict = tmp_obj_dict["seen"]
    data_frame = {}
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        if main_id.strip() == "":
            continue
        for tmp_row in data_frame["data"][main_id]:
            combo_id = main_id
            for f in combo_flist_one:
                val = tmp_row[tmp_fl.index(f)]
                val = val[:100] + " ... " if len(val) > 100 else val
                if in_file.find("proteoform_glycosylation_sites_") != -1 and f == "saccharide":
                    if val not in glycan_dict:
                        val = ""
                if f == "glycosylation_type":
                    val = val.lower()
                combo_id += "|" + val
            combo_id = combo_id.strip()
            if combo_id not in local_seen_dict:
                record_count += 1
                local_seen_dict[combo_id] = True
            
            obj_one = {}    
            for prop in prop_dict:
                f = prop_dict[prop]
                obj_one[prop] = tmp_row[tmp_fl.index(f)] if f in tmp_fl else ""
                if prop in ["do_id", "doid"]:
                    xref_key = tmp_row[tmp_fl.index(xref_info[0])]
                    xref_id = tmp_row[tmp_fl.index(xref_info[1])]
                    do_id = combo_id.split("|")[1]
                    obj_one[prop] = do_id
                    if do_id == "":
                        combo_id = "%s|%s-%s" % (main_id, xref_key.split("_")[-1], xref_id)
                        obj_one[prop] = "000000:%s-%s" % (xref_key.split("_")[-1], xref_id)
                        database_label = tmp_row[tmp_fl.index("database_label")]
                        do_name = "%s [%s disease name]" % (database_label, xref_badge)
                        obj_one["name"] = do_name
                        obj_one["url"] = libgly.get_xref_url(map_dict, "protein_xref_do_placeholder", "",is_cited)
                    else:
                        do_name = ""
                        if do_id in disease_id2names["do"]:
                            do_name = disease_id2names["do"][do_id]["name"]
                            do_name = do_name[0].upper() + do_name[1:]+" [DO disease name]"
                        obj_one["name"] = do_name
                        obj_one["url"] = libgly.get_xref_url(map_dict, "protein_xref_do", do_id,is_cited)

            if combo_id not in seen_dict:
                seen_dict[combo_id] = True
                tmp_obj_dict[combo_id] = obj_one
                if combo_flist_two != []:
                    obj_one["evidence"] = []


            if combo_flist_two != []:
                xref_key = tmp_row[tmp_fl.index(xref_info[0])]
                xref_id = tmp_row[tmp_fl.index(xref_info[1])]
                if xref_key in map_dict["xrefkey2badge"]:
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
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
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


def load_sequence_info(seq_dict, header_dict):

    file_list = glob.glob(data_dir + "/reviewed/*_protein_allsequences.fasta")
    for in_file in file_list:
        for record in SeqIO.parse(in_file, "fasta"):
            seq_id = record.id.split("|")[1]
            seq_dict[seq_id] = str(record.seq.upper())
   
     
    file_list = glob.glob(data_dir + "/reviewed/*_protein_sequenceinfo.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            isoform = row[f_list.index("uniprotkb_isoform_ac")]
            header_dict[isoform] = row[f_list.index("sequence_header")]

    return



def load_isoformlocus_info(isoformlocus_dict, species_name):

    ensembl_xref_key = custom_xref_key_dict[species_name] if species_name in custom_xref_key_dict else "protein_xref_ensembl"   
 
    s = species_list[0] if len(species_list) == 1 else "*"
    file_list = glob.glob(data_dir + "/reviewed/%s_protein_transcriptlocus.csv" % (s))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        tmp_fl = data_frame["fields"]
        for tmp_row in data_frame["data"]:
            isoform = tmp_row[tmp_fl.index("uniprotkb_isoform_ac")]
            transcript_id = tmp_row[tmp_fl.index("transcript_id")]
            peptide_id = tmp_row[tmp_fl.index("peptide_id")]
            t_url = libgly.get_xref_url(map_dict, ensembl_xref_key, transcript_id,is_cited)
            p_url = libgly.get_xref_url(map_dict, ensembl_xref_key, peptide_id,is_cited)
            start_pos = int(tmp_row[tmp_fl.index("start_pos")])
            end_pos = int(tmp_row[tmp_fl.index("end_pos")])
            strand = tmp_row[tmp_fl.index("strand")]
            strand_sign = "+"
            if strand == "0":
                start_pos = int(tmp_row[tmp_fl.index("end_pos")])
                end_pos = int(tmp_row[tmp_fl.index("start_pos")])
                strand_sign = "-"
            chr_id = tmp_row[tmp_fl.index("chromosome_id")]
            if isoform in isoformlocus_dict:
                continue
            isoformlocus_dict[isoform] = {
                "chromosome":chr_id
                ,"start_pos":start_pos
                ,"end_pos":end_pos
                ,"strand":strand_sign
                ,"evidence":[
                    {
                        "database":"Ensembl Transcript"
                        ,"id":transcript_id
                        ,"url":t_url
                    }
                    ,{
                        "database":"Ensembl Peptide"
                        ,"id":peptide_id
                        ,"url":p_url
                    }
                ]
            }
    return


def load_synthesized_glycans(synthesized_glycans, glycan_class):

    seen_glytoucan = {}
    in_file = data_dir + "/reviewed/glycan_enzyme.csv"
    sheet_obj = {}
    libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "glytoucan_ac")
    tmp_fl = sheet_obj["fields"]
    for glytoucan_ac in sheet_obj["data"]:
        for row in sheet_obj["data"][glytoucan_ac]:
            canon = row[tmp_fl.index("uniprotkb_canonical_ac")]
            if canon not in synthesized_glycans:
                synthesized_glycans[canon] = []
                seen_glytoucan[canon] = {}
            g_type, g_subtype = "Other", "Other"
            if glytoucan_ac in glycan_class:
                g_type = glycan_class[glytoucan_ac]["type"]
                g_subtype = glycan_class[glytoucan_ac]["subtype"]
            if glytoucan_ac not in seen_glytoucan[canon]:
                o = {"glytoucan_ac":glytoucan_ac, "type":g_type, "subtype":g_subtype}
                synthesized_glycans[canon].append(o)
                seen_glytoucan[canon][glytoucan_ac] = True
    return

def load_glycan_class(glycan_class):


    data_frame = {}
    in_file = data_dir + "/reviewed/glycan_classification.csv"
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        for tmp_row in data_frame["data"][main_id]:
            g_type = tmp_row[tmp_fl.index("glycan_type")].strip()
            g_subtype = tmp_row[tmp_fl.index("glycan_subtype")].strip()
            g_type = "Other" if g_type == "" else g_type
            g_subtype = "Other" if g_subtype == "" else g_subtype
            glycan_class[main_id] = {"type":g_type, "subtype":g_subtype}

    return 


def load_glycan_masterlist(glycan_dict):

    is_motif = {}
    in_file = data_dir +  "/reviewed/glycan_motif.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("motif_ac_xref")]
        is_motif[ac] = True

    data_frame = {}
    in_file = data_dir +  "/reviewed/glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        glycan_dict[ac] = True


    return 


def load_genelocus_info(genelocus_dict, species_name):

    ensembl_xref_key = custom_xref_key_dict[species_name] if species_name in custom_xref_key_dict else "protein_xref_ensembl"
    file_list = glob.glob(data_dir + "/reviewed/%s_protein_genelocus.csv" % (species_name))                
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        tmp_fl = data_frame["fields"]
        for tmp_row in data_frame["data"]:
            canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
            ensembl_gene_id = tmp_row[tmp_fl.index("ensembl_gene_id")]
            start_pos = int(tmp_row[tmp_fl.index("start_pos")])
            end_pos = int(tmp_row[tmp_fl.index("end_pos")])
            strand = tmp_row[tmp_fl.index("strand")]
            gene_symbol = tmp_row[tmp_fl.index("gene_symbol")]
            strand_sign = "+"
            if strand == "0":
                start_pos = int(tmp_row[tmp_fl.index("end_pos")])
                end_pos = int(tmp_row[tmp_fl.index("start_pos")])
                strand_sign = "-"
            chr_id = tmp_row[tmp_fl.index("chromosome_id")]
            gene_url = libgly.get_xref_url(map_dict, ensembl_xref_key, ensembl_gene_id,is_cited)
            if canon not in genelocus_dict:
                genelocus_dict[canon] = {
                    "chromosome":chr_id
                    ,"start_pos":start_pos
                    ,"end_pos":end_pos
                    ,"strand":strand_sign
                    ,"evidence":[
                        {
                            "database":"Ensembl Gene"
                            ,"id":ensembl_gene_id
                            ,"url":gene_url
                        }
                    ]
                }
    return




def load_cluster_info(cls_dict, recname_dict, submittedname_dict):

    canon2genename = {}
    for in_file in glob.glob(data_dir + "/reviewed/*_protein_masterlist.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        tmp_fl = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[tmp_fl.index("uniprotkb_canonical_ac")]
            canon2genename[canon] = row[tmp_fl.index("gene_name")]

    homolog_dict = {}
    member_dict = {}
    data_frame = {}
    in_file = data_dir + "/reviewed/protein_homolog_clusters.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    tmp_fl = data_frame["fields"]
    for row in data_frame["data"]:
        homolog_cluster_id = row[tmp_fl.index("homolog_cluster_id")]
        canon = row[tmp_fl.index("uniprotkb_canonical_ac")]
        if canon not in seq_dict:
            continue
        xref_key = row[tmp_fl.index("xref_key")]
        xref_id = row[tmp_fl.index("xref_id")]
        cls_id = "%s %s %s" % (homolog_cluster_id, xref_key, xref_id)
        recname = ""
        if canon in recname_dict["protein"]:
            recname =  recname_dict["protein"][canon]
        elif canon in submittedname_dict["protein"]:
            recname = submittedname_dict["protein"][canon]
        tax_id = row[tmp_fl.index("tax_id")]
        tax_name = species_obj[str(tax_id)]["long_name"]
        common_name = species_obj[str(tax_id)]["common_name"]
        glygen_name = species_obj[str(tax_id)]["glygen_name"]
        gene_name = canon2genename[canon] if canon in canon2genename else ""
        if canon not in homolog_dict:
            homolog_dict[canon] = {
                "uniprot_canonical_ac":canon
                ,"protein_name":recname
                ,"gene_name":gene_name
                ,"tax_id":int(tax_id)
                ,"organism":tax_name
                ,"sequence": {
                    "sequence": seq_dict[canon]
                    ,"length": len(seq_dict[canon])
                }
                ,"evidence":[]
            }
            for name_field in ["glygen_name", "common_name"]:
                name_value = species_obj[str(tax_id)][name_field]
                if name_value != "":
                    homolog_dict[canon][name_field] = name_value

        if cls_id not in member_dict:
            member_dict[cls_id] = {}
        member_dict[cls_id][canon] = True

    evdn_dict = {}
    for cls_id in member_dict:
        member_list = member_dict[cls_id]
        for canon_one in member_list:
            for canon_two in member_list:
                if canon_one == canon_two:
                    continue
                if canon_one not in evdn_dict:
                    evdn_dict[canon_one] = {}
                if canon_two not in evdn_dict[canon_one]:
                    evdn_dict[canon_one][canon_two] = {}
                if cls_id not in evdn_dict[canon_one][canon_two]:
                    evdn_dict[canon_one][canon_two][cls_id] = True

    for canon_one in evdn_dict:
        for canon_two in evdn_dict[canon_one]:
            combo_id = "%s|%s" % (canon_one, canon_two)
            if combo_id not in main_dict["orthologs"]:
                main_dict["orthologs"][combo_id] = json.loads(json.dumps(homolog_dict[canon_two]))
            for cls_id in evdn_dict[canon_one][canon_two].keys():
                xref_key, xref_id = cls_id.split(" ")[1], cls_id.split(" ")[2]
                if xref_key == "protein_xref_oma":
                    xref_id = canon_one.split("-")[0]
                xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                o = {"id":xref_id, "url":xref_url, "database":xref_badge}
                main_dict["orthologs"][combo_id]["evidence"].append(o)


    return


def load_name_sheet_info():

    return {
        "protein_recnames":{
            "badge":"UniProtKB",
            "resource":"uniprotkb", 
            "xref_key":"protein_xref_uniprotkb_proteinname",
            "fieldlist":["recommended_name_full","recommended_name_short","ec_name"]
        },
        "protein_altnames":{
            "badge":"UniProtKB",
            "resource":"uniprotkb",
            "xref_key":"protein_xref_uniprotkb_proteinname",
            "fieldlist":["alternative_name_full","alternative_name_short","ec_name"]
        },
        "protein_submittednames":{
            "badge":"UniProtKB",
            "resource":"uniprotkb",
            "xref_key":"protein_xref_uniprotkb_proteinname",
            "fieldlist":["submitted_name_full","submitted_name_short","ec_name"]
        },
        "protein_proteinnames_refseq":{
            "badge":"RefSeq",
            "resource":"refseq",
            "xref_key":"protein_xref_refseq",
            "fieldlist":["refseq_protein_name"]
        },
        "protein_genenames_uniprotkb":{
            "badge":"UniProtKB",
            "resource":"uniprotkb",
            "xref_key":"protein_xref_uniprotkb",
            "fieldlist":["gene_symbol_recommended","gene_symbol_alternative","orf_name"]
        },  
        "protein_genenames_refseq":{
            "badge":"RefSeq",
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


def get_disease_synonyms(do_id, mondo_id, mim_id):
    syn_obj_list = {
        "do":{"synlist":[], "xrefid":""},
        "mondo":{"synlist":[], "xrefid":""},
        "omim":{"synlist":[], "xrefid":""}
    }
    
    if do_id != "" and do_id in disease_id2names["do"]:
        syn_obj_list["do"]["xrefid"] = do_id
        syn_obj_list["do"]["xrefidlabel"] = "DOID:" + do_id
        syn_obj_list["do"]["synlist"] = [disease_id2names["do"][do_id]["name"]]
        syn_obj_list["do"]["synlist"] += disease_id2names["do"][do_id]["synonyms"]
        
    if mondo_id != "" and mondo_id in disease_id2names["mondo"]:
        syn_obj_list["mondo"]["xrefid"] = mondo_id
        syn_obj_list["mondo"]["xrefidlabel"] = "MONDO:" + mondo_id
        syn_obj_list["mondo"]["synlist"] = [disease_id2names["mondo"][mondo_id]["name"]]
        syn_obj_list["mondo"]["synlist"] += disease_id2names["mondo"][mondo_id]["synonyms"]
        
    if  mim_id != "" and mim_id in disease_id2names["omim"]:
        syn_obj_list["omim"]["xrefid"] = mim_id
        syn_obj_list["omim"]["xrefidlabel"] = "MIM:" + mim_id
        syn_obj_list["omim"]["synlist"] = [disease_id2names["omim"][mim_id]["name"]]
        syn_obj_list["omim"]["synlist"] += disease_id2names["omim"][mim_id]["synonyms"]

    return syn_obj_list


def get_disease_info(do_id, mondo_id, mim_id):
    disease_info = {}
    json_file = ""
    if do_id != "":
        json_file = "jsondb/diseasedb/doid.%s.json" % (do_id.replace("DOID:", ""))
    elif mondo_id != "":
        json_file = "jsondb/diseasedb/mondo.%s.json" % (mondo_id.replace("MONDO:", ""))
    elif mim_id != "":
        json_file = "jsondb/diseasedb/mim.%s.json" % (mim_id.replace("OMIM:", ""))

    if os.path.isfile(json_file):
        disease_info = json.loads(open(json_file, "r").read())


    return disease_info



def get_disease_info_old(do_id, mondo_id, mim_id):

    disease_info = {}
    if do_id != "":
        name_xref_key = "protein_xref_do"
        name_xref_id = do_id
        name_xref_url = libgly.get_xref_url(map_dict, name_xref_key, name_xref_id,is_cited)
        name_xref_badge = libgly.get_xref_badge(map_dict, name_xref_key)
        name = disease_id2names["do"][name_xref_id]["name"] if name_xref_id in disease_id2names["do"] else ""
        disease_desc = disease_id2desc["do"][name_xref_id] if name_xref_id in disease_id2desc["do"] else ""
        disease_info = {
                "key":name_xref_key,"id":"DOID:"+name_xref_id, "badge":name_xref_badge, 
            "url":name_xref_url, "name":name, "description":disease_desc
        }
    elif mondo_id != "":
        name_xref_key = "protein_xref_mondo"
        name_xref_id = mondo_id
        name_xref_url = libgly.get_xref_url(map_dict, name_xref_key, name_xref_id,is_cited)
        name_xref_badge = libgly.get_xref_badge(map_dict, name_xref_key)
        name = disease_id2names["mondo"][name_xref_id]["name"]
        disease_desc = disease_id2desc["mondo"][name_xref_id]
        disease_info = {
                "key":name_xref_key,"id":"MONDO:"+name_xref_id, "badge":name_xref_badge,
            "url":name_xref_url, "name":name, "description":disease_desc
        }
    elif mim_id != "" and mim_id in disease_id2names["omim"]:
        name_xref_key = "protein_xref_omim"
        name_xref_id = mim_id
        name_xref_url = libgly.get_xref_url(map_dict, name_xref_key, name_xref_id,is_cited)
        name_xref_badge = libgly.get_xref_badge(map_dict, name_xref_key)
        name = disease_id2names["omim"][name_xref_id]["name"]
        disease_desc = disease_id2desc["omim"][name_xref_id]
        disease_info = {
                "key":name_xref_key,"id":"MIM:"+name_xref_id, "badge":name_xref_badge, 
            "url":name_xref_url, "name":name, "description":disease_desc
        }

    disease_obj = {}
    if disease_info != {}:
        disease_name = disease_info["name"]
        disease_desc = disease_info["description"]

        name_xref_key, name_xref_id = disease_info["key"], disease_info["id"]
        name_xref_badge, name_xref_url = disease_info["badge"],disease_info["url"]
        disease_synonyms = []
        syn_obj_list = get_disease_synonyms(do_id, mondo_id, mim_id)
        seen_syn = {}
        for k in syn_obj_list:
            syn_xref_key = "protein_xref_" + k
            syn_xref_id = syn_obj_list[k]["xrefid"]
            syn_xref_url = libgly.get_xref_url(map_dict,syn_xref_key,syn_xref_id,is_cited)
            syn_xref_badge = libgly.get_xref_badge(map_dict, syn_xref_key)
            for syn in syn_obj_list[k]["synlist"]:
                if syn.lower() == disease_name.lower():
                    continue
                syn_combo_id = "%s|%s|%s" % (syn_xref_key,syn_xref_id,syn.lower())
                if syn_combo_id not in seen_syn:
                    o = {"name":syn, "description":"", "id": syn_obj_list[k]["xrefidlabel"],
                                        "resource":syn_xref_badge,"url": syn_xref_url}
                    disease_synonyms.append(o)
                    seen_syn[syn_combo_id] = True
        disease_obj = {
            "recommended_name":{
                "name":disease_name, "description":disease_desc, "id": name_xref_id,
                "resource":name_xref_badge,"url": name_xref_url
            }
            ,"synonyms":disease_synonyms
        }

    return disease_obj

def load_refseq_dict():

    refseq_dict = {}
    for combo_id in main_dict["protein_names"]:
        canon = combo_id.split("|")[0]
        if combo_id.split("|")[-1] == "refseq":
            obj = main_dict["protein_names"][combo_id]
            if canon not in refseq_dict:
                ac = obj["id"]
                url = libgly.get_xref_url(map_dict, "protein_xref_refseq", ac,is_cited)
                refseq_dict[canon] = {"ac":ac, "name":obj["name"], "url":url}
    return refseq_dict

def load_is_reported(expected_dslist):


    #predicted_xref_key_list = [ "protein_xref_uniprotkb_gly"]
    predicted_list = [
        "proteoform_glycosylation_sites_uniprotkb", 
        "proteoform_glycosylation_sites_predicted_isoglyp"
    ]    

    is_reported = {}
    for file_name in expected_dslist:
        file_ext = "fasta" if file_name.find("protein_allsequences") != -1 else "csv"
        in_file = "%s/reviewed/%s.%s" % (data_dir,file_name,file_ext)
        predicted_flag = False
        for val in predicted_list:
            if file_name.find(val) != -1:
                predicted_flag = True
                break
        if file_name.find("proteoform_glycosylation_sites_") != -1:
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
                main_id = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                aa_pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                glytoucan_ac = tmp_row[tmp_fl.index("saccharide")]
                start_pos = tmp_row[tmp_fl.index("start_pos")]
                end_pos = tmp_row[tmp_fl.index("end_pos")]
                xref_key = tmp_row[tmp_fl.index("xref_key")]
                #if xref_key not in predicted_xref_key_list:
                if predicted_flag == False:
                    combo_id = "%s|%s|%s|%s" % (main_id,start_pos, end_pos,glytoucan_ac)
                    is_reported[combo_id] = True
                    combo_id = "%s|%s|%s|%s" % (main_id,start_pos, end_pos,"")
                    is_reported[combo_id] = True
            FR.close()

    return is_reported




def make_objects(species_name):


    # set order to 1000 to move datasets to the end
    # set order to 0 to move datasets to the start
    order_dict = {"_protein_masterlist":0, "_protein_expression_disease": 1000}
    tmp_dslist =  csvutil.get_expected_dslist("protein", [], order_dict)
    expected_dslist = []
    for f in tmp_dslist:
        sp, mol = f.split("_")[0], f.split("_")[1]
        if mol in ["protein", "proteoform"]:
            if sp == species_name:
                expected_dslist.append(f)
        else:
            expected_dslist.append(f)

    #print (json.dumps(expected_dslist, indent=4))
    #print (len(expected_dslist))
    #exit()

    missing_list = []
    for file_name in expected_dslist:
        in_file = "%s/reviewed/%s.csv" % (data_dir,file_name)
        if os.path.isfile(in_file) == False:
            missing_list.append(in_file)
    if missing_list != []:
        print ("Dataset files missing:")
        print (json.dumps(missing_list, indent=4))
        sys.exit()

    
    sec_info = json.loads(open("generated/misc/protein_sectioninfo.json", "r").read())
    for sec in sec_info:
        main_dict[sec] = {"seen":{}}


    log_file = "logs/make-proteindb.log"


    msg = "make-proteindb: loading dictionaries for %s" % (species_name)
    csvutil.write_log_msg(log_file, msg, "a")
    csvutil.load_dictionaries(map_dict, "generated/misc/")


    msg = "make-proteindb: loading protein names for %s" % (species_name)
    csvutil.write_log_msg(log_file, msg, "a")
    recname_dict = {"protein":{}, "gene":{}}
    submittedname_dict = {"protein":{}}
    load_protein_names(recname_dict, submittedname_dict)

    msg = "make-proteindb: loading enzyme_ann for %s" % (species_name)
    csvutil.write_log_msg(log_file, msg, "a")
    enzyme_ann = {}
    load_enzyme_ann(enzyme_ann, species_name)

    msg = "make-proteindb: loading enzyme_dict for %s" % (species_name)
    csvutil.write_log_msg(log_file, msg, "a")
    enzyme_dict = {}
    load_enzyme_dict(enzyme_dict, species_name)

    
    msg = "make-proteindb: loading gene locus"
    csvutil.write_log_msg(log_file, msg, "a")
    genelocus_dict = {}
    load_genelocus_info(genelocus_dict, species_name)
   
    msg = "make-proteindb: loading isoform info"
    csvutil.write_log_msg(log_file, msg, "a")
    isoformlocus_dict = {}
    load_isoformlocus_info(isoformlocus_dict, species_name)
 
    msg = "make-proteindb: loading reaction_dict participant_dict ac2rxn for %s" % (species_name)
    csvutil.write_log_msg(log_file, msg, "a")
    reaction_dict, participant_dict, ac2rxn = {}, {}, {}
    load_reactions(reaction_dict,participant_dict, ac2rxn, enzyme_dict, species_name)

    msg = "make-proteindb: loading pathway_dict"
    csvutil.write_log_msg(log_file, msg, "a")
    pathway_dict = {}
    load_pathways(pathway_dict)

    msg = "make-proteindb: loading disease_id2names, disease_id2desc"
    csvutil.write_log_msg(log_file, msg, "a")
    load_disease_names()


    msg = "make-proteindb: loading refseq info"
    csvutil.write_log_msg(log_file, msg, "a")
    refseq_dict = load_refseq_dict()

    msg = "make-proteindb: loading cls_dict"
    csvutil.write_log_msg(log_file, msg, "a")
    cls_dict = {}
    load_cluster_info(cls_dict, recname_dict, submittedname_dict)

    msg = "make-proteindb: loading is_reported"
    csvutil.write_log_msg(log_file, msg, "a")
    is_reported = load_is_reported(expected_dslist)


    biomarker_dict = csvutil.get_biomarker_dict("protein")

    seen_combo_glycosylation = {}
    seen_combo_disease = {}

    file_idx = 1
    file_count = len(expected_dslist)
    for file_name in expected_dslist:
        file_ext = "fasta" if file_name.find("protein_allsequences") != -1 else "csv" 
        in_file = "%s/reviewed/%s.%s" % (data_dir,file_name,file_ext)
        if os.path.isfile(in_file) == False:
            msg = "make-proteindb: file %s does NOT exist!" % (in_file)
            csvutil.write_log_msg(log_file, msg, "a")
            sys.exit()
  
        msg = "make-proteindb: %s.%s [%s/%s]" % (file_name,file_ext, file_idx, file_count)
        csvutil.write_log_msg(log_file, msg, "a")
        file_idx += 1
        
        #--> uniprot_canonical_ac, uniprot_ac, species
        sheet_name = "protein_masterlist"
        if file_name.find(sheet_name) != -1:
            species = file_name.split("_")[0]
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
            tmp_fl = data_frame["fields"]
            for main_id in data_frame["data"]:
                if main_id in main_dict["uniprot_canonical_ac"]:
                    continue
                uniprotkb_ac = main_id.split("-")[0]
               
                prop_name = "uniprot_canonical_ac"
                main_dict[prop_name][main_id] = main_id
                combo_list = ["uniprotkb_canonical_ac"]
                update_record_stat(record_stat,  file_name, prop_name, 1,combo_list)
                
                prop_name = "uniprot_ac"
                main_dict[prop_name][main_id] = uniprotkb_ac
                combo_list = ["uniprotkb_canonical_ac", "uniprot_ac"]
                update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)


                prop_name = "reactions"
                participant_list, enzyme_list = [], []
                if main_id in ac2rxn:
                    for rxn_id in ac2rxn[main_id]:
                        combo_id = "%s|%s" % (main_id, rxn_id)
                        main_dict[prop_name][combo_id]= reaction_dict[rxn_id]
                        for k in ["output_participant_list", "input_participant_list"]:
                            if k in reaction_dict[rxn_id]:
                                for p_id in reaction_dict[rxn_id][k]:
                                    participant_list.append(p_id)
                        k = "enzyme_list"
                        if k in reaction_dict[rxn_id]:
                            for e_id in reaction_dict[rxn_id][k]:
                                if e_id in enzyme_dict:
                                    enzyme_list.append(e_id)

                prop_name = "reaction_participants"
                for p_id in participant_list:
                    combo_id = "%s|%s" % (main_id, p_id)
                    main_dict[prop_name][combo_id]= participant_dict[p_id]
                
                prop_name = "reaction_enzymes"
                for e_id in enzyme_list:
                    combo_id = "%s|%s" % (main_id, e_id)
                    main_dict[prop_name][combo_id] = enzyme_dict[e_id]




                seq_header = header_dict[main_id] if main_id in header_dict else ""
                o = {
                    "sequence":seq_dict[main_id]
                    ,"length":len(seq_dict[main_id])
                    ,"header":seq_header
                }
                prop_name = "sequence"
                main_dict[prop_name][main_id] = o
                combo_list = ["uniprotkb_canonical_ac", "sequence"]
                update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)

                xref_key, xref_id = "protein_xref_uniprotkb", uniprotkb_ac
                gene_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                tmp_obj = {}
                locus_obj = genelocus_dict[main_id] if main_id in genelocus_dict else tmp_obj
                tmp_row = data_frame["data"][main_id][0]
                gene_name = tmp_row[tmp_fl.index("gene_name")]
                o = { "name":gene_name ,"url":gene_url,"locus":locus_obj}
                combo_id = "%s|%s" % (main_id, gene_name)
                prop_name = "gene"
                main_dict[prop_name][combo_id] = o
                combo_list = ["uniprotkb_canonical_ac","gene_symbol"]
                update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)

                
                url = libgly.get_xref_url(map_dict, "protein_xref_uniprotkb", uniprotkb_ac,is_cited)
                o = {
                    "name":species_obj[species]["long_name"]
                    ,"taxid":species_obj[species]["tax_id"]
                    ,"evidence":[
                        {"database":"UniProtKB", "id":uniprotkb_ac,"url":url}
                    ]   
                }
                for nm in ["common_name", "glygen_name"]:
                    name_val = species_obj[str(species_obj[species]["tax_id"])][nm]
                    if name_val != "":
                        o[nm] = name_val
                o["reference_species"] = "%s [%s]" % (o["name"], o["taxid"])

                combo_id = "%s|%s" % (main_id, species_obj[species]["tax_id"])
                prop_name = "species"
                main_dict[prop_name][combo_id] = o
                combo_list = ["uniprotkb_canonical_ac","tax_id"]
                update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)
                
                isoform_list = []
                for tmp_row in data_frame["data"][main_id]:
                    for f in ["reviewed_isoforms", "unreviewed_isoforms"]:
                        if tmp_row[tmp_fl.index(f)] != "":
                            isoform_list.append(tmp_row[tmp_fl.index(f)])
                for isoform in isoform_list:
                    locus_obj = {}
                    if isoform in isoformlocus_dict:
                        locus_obj = isoformlocus_dict[isoform]
                    xref_key = "protein_xref_uniprot_isoform"
                    isoform_url = ""
                    if xref_key in map_dict["xrefkey2url"]:
                        isoform_url = map_dict["xrefkey2url"][xref_key][0] % (uniprotkb_ac,isoform)
                    isoform_seq, isoform_header = "", ""
                    if isoform in seq_dict:
                        isoform_seq, isoform_header = seq_dict[isoform], header_dict[isoform]
                    o = {
                        "isoform_ac":isoform
                        ,"url":isoform_url
                        ,"sequence":{
                            "sequence":isoform_seq
                            ,"length":len(isoform_seq)
                            ,"header":isoform_header
                        },
                        "locus":locus_obj
                    }
                    combo_id = "%s|%s" % (main_id, isoform)
                    prop_name = "isoforms"
                    main_dict[prop_name][combo_id] = o
                    combo_list = ["uniprotkb_canonical_ac","isoforms"]
                    update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)


        #--> mass, uniport_id
        sheet_name = "protein_info_uniprotkb"
        if file_name.find(sheet_name) != -1:
            species = file_name.split("_")[0]
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
            tmp_fl = data_frame["fields"]
            for main_id in data_frame["data"]:
                for tmp_row in data_frame["data"][main_id][:1]:
                    protein_mass = float(tmp_row[tmp_fl.index("uniprotkb_protein_mass")])
                    protein_length = int(tmp_row[tmp_fl.index("uniprotkb_protein_length")])
                    uniprotkb_id = tmp_row[tmp_fl.index("uniprotkb_id")]
                    prop_name = "uniprot_id"
                    main_dict[prop_name][main_id] = uniprotkb_id
                    combo_list = ["uniprotkb_canonical_ac","uniprotkb_id"]
                    update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)
                    
                    prop_name = "mass"
                    main_dict[prop_name][main_id] = {"chemical_mass": protein_mass, 
                            "monoisotopic_mass":protein_mass}
                    combo_list = ["uniprotkb_canonical_ac","uniprotkb_protein_mass"]
                    update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)

        #--> refseq
        sheet_name = "protein_info_refseq"
        if file_name.find(sheet_name) != -1:
            species = file_name.split("_")[0]
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
            tmp_fl = data_frame["fields"]
            for main_id in data_frame["data"]:
                for tmp_row in data_frame["data"][main_id][:1]:
                    refseq_ac = tmp_row[tmp_fl.index("p_refseq_ac_best_match")]
                    refseq_proteinname = tmp_row[tmp_fl.index("refseq_protein_name")]
                    refseq_url = libgly.get_xref_url(map_dict,"protein_xref_refseq", refseq_ac,is_cited)
                    refseq_summary = tmp_row[tmp_fl.index("refseq_protein_summary")]
                    prop_name = "refseq"
                    main_dict[prop_name][main_id] = {
                            "ac":refseq_ac,"name":refseq_proteinname, 
                            "summary":refseq_summary, "url": refseq_url 
                    }
                    combo_list = ["uniprotkb_canonical_ac","refseq"]
                    update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)



             

        #--> keywords
        for sheet_name in ["protein_glycosyltransferase", "protein_glycohydrolase"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    if main_id not in main_dict["keywords"]:
                        main_dict["keywords"][main_id] = []
                    kw_list = ["enzyme", sheet_name.split("_")[-1] + "-activity"]
                    for kw in kw_list:
                        if kw not in main_dict["keywords"][main_id]:
                            prop_name = "keywords"
                            main_dict[prop_name][main_id].append(kw)
                            combo_list = ["uniprotkb_canonical_ac","keywords"]
                            update_record_stat(record_stat,  file_name, prop_name, 1, combo_list)

        #--> interactions
        for sheet_name in ["protein_matrixdb"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "interactions"
                prop_dict = {"interactor_id":"saccharide", "interaction_name":"matrix_db_label",
                    "interaction_type":"interaction_type"
                }
                xref_info = ["xref_key", "xref_id"]
                combo_flist_one = ["saccharide"]
                combo_flist_two = ["xref_key", "xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj,in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    if combo_id.split("|")[1] == "":
                        load_obj[combo_id]["deleteflag"] = True

        #--> publication
        for sheet_name in ["protein_citations_", "proteoform_citations_"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "publication"
                prop_dict = {"title":"title", 
                        "journal":"journal_name","date":"publication_date","authors":"authors"}
                xref_info = ["src_xref_key", "src_xref_id"]
                combo_flist_one = ["xref_key", "xref_id"]
                combo_flist_two = ["src_xref_key", "src_xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj,in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id, xref_key, xref_id = combo_id.split("|")
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


        #--> structures
        for sheet_name in ["protein_pdb_map"]:
            if file_name.find(sheet_name) != -1 and file_name.split("_")[0] != "protein":
                species = file_name.split("_")[0]
                prop_name = "structures"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        sequence_region = tmp_row[tmp_fl.index("sequence_region")]
                        pdb_id = tmp_row[tmp_fl.index("pdb_id")]
                        pdb_chain = tmp_row[tmp_fl.index("pdb_chain")]
                        start_pos = tmp_row[tmp_fl.index("start_pos")]
                        end_pos = tmp_row[tmp_fl.index("end_pos")]
                        overlap = tmp_row[tmp_fl.index("overlap_ratio")]
                        method = tmp_row[tmp_fl.index("experimental_method")]
                        resolution = tmp_row[tmp_fl.index("resolution")]
                        selection_flag = tmp_row[tmp_fl.index("selection_flag")]
                        start_pos = int(start_pos) if start_pos.isdigit() else -1
                        end_pos = int(end_pos) if end_pos.isdigit() else -1
                        resolution = float(resolution) if resolution.isdigit() else -1.0
                        pdb_type = "alphafold" if method in ["AlphaFold"] else "experimental"
                        if selection_flag == "True" and pdb_id in pdbid2url:
                            combo_id_one = "%s|%s" % (main_id,sequence_region)
                            load_obj[combo_id_one] = {
                                "pdb_id":pdb_id,
                                "pdb_chain":pdb_chain,
                                "start_pos":start_pos,
                                "end_pos":end_pos,
                                "overlap":overlap,
                                "method":method,
                                "type":pdb_type,
                                "resolution":resolution,
                                "url":pdbid2url[pdb_id]                    
                            }        


 
        #--> disease
        for sheet_name in ["protein_disease"]:
            if file_name.find(sheet_name) != -1 and file_name.split("_")[0] != "protein":
                species = file_name.split("_")[0]
                prop_name = "disease"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        xref_key = tmp_row[tmp_fl.index("xref_key")]
                        xref_id = tmp_row[tmp_fl.index("xref_id")]
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                        xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                        do_id = tmp_row[tmp_fl.index("do_id")].strip()
                        mondo_id = tmp_row[tmp_fl.index("mondo_id")].strip()
                        mim_id = xref_id if xref_key == "protein_xref_omim" else ""
                        name_xref_id, name_xref_url, name_xref_badge = "", "", ""
                    
                        disease_obj = get_disease_info(do_id, mondo_id, mim_id)
                        if disease_obj == {}:
                            continue
                        name_xref_id = disease_obj["recommended_name"]["id"]
                        name_xref_badge = disease_obj["recommended_name"]["resource"]
                        combo_id_one = "%s|%s|%s" % (main_id,name_xref_badge,name_xref_id)
                        if combo_id_one not in seen_combo_disease:
                            load_obj[combo_id_one] = disease_obj
                            load_obj[combo_id_one]["evidence"] = []
                            seen_combo_disease[combo_id_one] = []
                        ev_obj = {"id":xref_id, "database":xref_badge, "url":xref_url}
                        combo_id_two = "%s|%s" % (xref_key,xref_id)
                        if combo_id_two not in seen_combo_disease[combo_id_one]:
                            load_obj[combo_id_one]["evidence"].append(ev_obj)
                            seen_combo_disease[combo_id_one].append(combo_id_two)






        #--> function
        for sheet_name in ["protein_function_uniprotkb", "protein_function_refseq"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "function"
                load_obj = main_dict[prop_name]
                prop_dict = {"annotation":"annotation"}
                xref_info = ["xref_key", "xref_id"]
                #combo_flist_one = ["annotation", "xref_key", "xref_id"]
                combo_flist_one = ["annotation"]
                combo_flist_two = ["xref_key", "xref_id"]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)

        #--> phosphorylation
        for sheet_name in ["proteoform_phosphorylation_sites_"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "phosphorylation"
                prop_dict = {
                    "start_pos":"phosphorylation_site_uniprotkb",
                    "end_pos":"phosphorylation_site_uniprotkb",
                    "kinase_uniprot_canonical_ac":"kinase_uniprotkb_canonical_ac",
                    "kinase_gene_name":"kinase_gene_name",
                    "residue":"amino_acid",
                    "comment":"phosphorylation_annotation"
                }
                xref_info = ["xref_key", "xref_id"]
                combo_flist_one = ["phosphorylation_site_uniprotkb",
                    "kinase_uniprotkb_canonical_ac", "amino_acid"
                ]
                combo_flist_two = ["xref_key", "xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)

                for combo_id in load_obj:
                    if "start_pos" in load_obj[combo_id]:
                        aa_pos = int(load_obj[combo_id]["start_pos"])
                        residue = load_obj[combo_id]["residue"]
                        load_obj[combo_id]["start_pos"]  = int(load_obj[combo_id]["start_pos"])
                        load_obj[combo_id]["end_pos"]  = int(load_obj[combo_id]["end_pos"])
                        load_obj[combo_id]["site_lbl"] =  "%s%s" % (residue, aa_pos)



         
        for sheet_name in ["proteoform_glycosylation_sites_"]:
            if file_name.find(sheet_name) != -1:
                prd_tool = ""
                if file_name.find("glycosylation_sites_predicted_") != -1:
                    prd_tool = file_name.split("_predicted_")[-1].replace(".csv", "")
                
                species = file_name.split("_")[0]
                prop_name = "glycosylation"
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
                        main_id = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                        tmp_aa_pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                        tmp_glytoucan_ac = tmp_row[tmp_fl.index("saccharide")]
                        tmp_amino_acid = tmp_row[tmp_fl.index("amino_acid")]
                        site_seq = tmp_row[tmp_fl.index("site_seq")] if "site_seq" in tmp_fl else ""
                        start_pos = tmp_row[tmp_fl.index("start_pos")]
                        end_pos = tmp_row[tmp_fl.index("end_pos")]
                        start_aa = tmp_row[tmp_fl.index("start_aa")]
                        end_aa = tmp_row[tmp_fl.index("end_aa")]
                        mining_tool = tmp_row[tmp_fl.index("mining_tools")] if "mining_tools" in tmp_fl else ""
                        aa_pos_list = tmp_aa_pos.split("|")
                        glytoucan_ac_list = tmp_glytoucan_ac.split("|")
                        amino_acid_list = tmp_amino_acid.split("|")
                        aa_pos = aa_pos_list[0]
                        glytoucan_ac = glytoucan_ac_list[0]
                        amino_acid = amino_acid_list[0]
                        if glytoucan_ac != "" and glytoucan_ac not in glycan_dict:
                            continue

                        gly_type = tmp_row[tmp_fl.index("glycosylation_type")]
                        if gly_type != "":
                            gly_type = gly_type[0].upper() + gly_type[1:]
                        gly_subtype = ""
                        if "glycosylation_subtype" in tmp_fl:
                            gly_subtype = tmp_row[tmp_fl.index("glycosylation_subtype")]
                        
                        xref_key = tmp_row[tmp_fl.index("xref_key")]
                        xref_id = tmp_row[tmp_fl.index("xref_id")]
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                        xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                        
                        src_xref_key = tmp_row[tmp_fl.index("src_xref_key")]
                        src_xref_id = tmp_row[tmp_fl.index("src_xref_id")]
                        src_xref_url = libgly.get_xref_url(map_dict, src_xref_key, src_xref_id,is_cited)
                        src_xref_badge = libgly.get_xref_badge(map_dict, src_xref_key)
                        
                        combo_id = "%s|%s|%s|%s" % (main_id,start_pos, end_pos,glytoucan_ac)
                        site_cat_one, site_cat_two = "", ""
                        if combo_id in is_reported:
                            site_cat_one = "reported"
                            site_cat_two = "reported"
                            if mining_tool != "":
                                site_cat_two = "automatic_literature_mining"
                            if glytoucan_ac != "":
                                site_cat_one = "reported_with_glycan"
                                site_cat_two += "_with_glycan"
                        else:
                            site_cat_one = "predicted"
                            site_cat_two = "predicted"
                            if glytoucan_ac != "":
                                site_cat_one = "predicted_with_glycan"
                                site_cat_two = "predicted_with_glycan"
                         

                        cond_list = []
                        #cond_list.append(aa_pos.strip() == "")
                        #cond_list.append(amino_acid == "")
                        cond_list.append(gly_type == "")
                        if True in cond_list:
                            continue
                        combo_id_one = "%s|%s|%s|%s|%s" % (main_id,start_pos,end_pos,glytoucan_ac,site_cat_one)
                        note_fields = ["uniprotkb_glycosylation_annotation_comment",
                            "glycosylation_subtype","curation_notes","additonal_notes"
                        ]
                        comment = ""
                        for f in note_fields:
                            comment += " " + tmp_row[tmp_fl.index(f)] if f in tmp_fl else ""
                        comment =  comment.strip()

                        prd_tool_list = []
                        ooo = {"label":prd_tool}
                        if prd_tool in tool_url_dict:
                            ooo["url"] = tool_url_dict[prd_tool]                  
                        prd_tool_list.append(ooo)

                        if combo_id_one not in seen_combo_glycosylation:
                            gly_obj = {
                                "glytoucan_ac":glytoucan_ac,
                                "type":gly_type,
                                "subtype":gly_subtype,
                                "prediction_tool":prd_tool,
                                "prediction_tool_list":prd_tool_list,
                                "site_category":site_cat_two,
                                "site_category_dict":{},
                                "site_seq":site_seq,
                                "relation":"attached",
                                "comment":comment
                            }
                            if start_pos.strip() != "":
                                gly_obj["start_pos"],gly_obj["start_aa"] = int(start_pos), start_aa
                            if end_pos.strip() != "":
                                gly_obj["end_pos"], gly_obj["end_aa"] = int(end_pos), end_aa
                                
                            if aa_pos.strip() != "" and amino_acid.strip() != "":
                                gly_obj["residue"] = amino_acid[0].upper()
                                if len(amino_acid) > 1:
                                    gly_obj["residue"] += amino_acid[1:].lower() 
                                gly_obj["site_lbl"] =  "%s%s" % (gly_obj["residue"], aa_pos)

                            if len(aa_pos_list) > 1:
                                gly_obj["alternate_start_pos_list"] = []
                                gly_obj["alternate_end_pos_list"] = []
                                for tmp_pos in aa_pos_list[1:]:
                                    gly_obj["alternate_start_pos_list"].append(int(tmp_pos))
                                    gly_obj["alternate_end_pos_list"].append(int(tmp_pos))
                            if len(glytoucan_ac_list) > 1:
                                gly_obj["alternate_glytoucan_ac_list"] = glytoucan_ac_list[1:]
                            if len(amino_acid_list) > 1:
                                gly_obj["alternate_amino_acid_list"] = amino_acid_list[1:]


                            load_obj[combo_id_one] = gly_obj
                            load_obj[combo_id_one]["evidence"] = []
                            seen_combo_glycosylation[combo_id_one] = []
                        
                        #overwrite conditions
                        cat_list_one = ["reported", "automatic_literature_mining", "predicted_with_glycan", "predicted"]
                        cat_list_two = ["automatic_literature_mining", "predicted_with_glycan", "predicted"]
                        cat_list_three = ["predicted_with_glycan", "predicted"]
                        cat_list_four = ["predicted"]
                        if site_cat_two == "reported_with_glycan" and load_obj[combo_id_one]["site_category"] in cat_list_one:
                            load_obj[combo_id_one]["site_category"] = site_cat_two
                        elif site_cat_two == "reported" and load_obj[combo_id_one]["site_category"] in cat_list_two:
                            load_obj[combo_id_one]["site_category"] = site_cat_two
                        elif site_cat_two == "automatic_literature_mining" and load_obj[combo_id_one]["site_category"] in cat_list_three:
                            load_obj[combo_id_one]["site_category"] = site_cat_two 
                        elif site_cat_two == "predicted_with_glycan" and load_obj[combo_id_one]["site_category"] in cat_list_four:
                            load_obj[combo_id_one]["site_category"] = site_cat_two
                           
                        load_obj[combo_id_one]["site_category_dict"][site_cat_two] = True
                        if mining_tool != "":
                            load_obj[combo_id_one]["mining_tool"] = mining_tool
                            load_obj[combo_id_one]["mining_tool_list"] = []
                            for tool in mining_tool.split(";"):
                                ooo = {"label":tool}
                                if tool in tool_url_dict:
                                    ooo["url"] = tool_url_dict[tool] 
                                load_obj[combo_id_one]["mining_tool_list"].append(ooo)


                        ev_obj = {"id":xref_id, "database":xref_badge, "url":xref_url}
                        combo_id_two = "%s|%s" % (xref_key,xref_id)
                        if combo_id_two not in seen_combo_glycosylation[combo_id_one]:
                            load_obj[combo_id_one]["evidence"].append(ev_obj)
                            seen_combo_glycosylation[combo_id_one].append(combo_id_two)
                        #xxxxxxx
                        ev_obj = {"id":src_xref_id, "database":src_xref_badge, "url":src_xref_url}
                        combo_id_two = "%s|%s" % (src_xref_key,src_xref_id)
                        if combo_id_two not in seen_combo_glycosylation[combo_id_one]:
                            load_obj[combo_id_one]["evidence"].append(ev_obj)
                            seen_combo_glycosylation[combo_id_one].append(combo_id_two)
                        
                        if main_id not in main_dict["keywords"]:
                            main_dict["keywords"][main_id] = []
                        kw_list = ["glycoprotein", "glycoprotein_" + site_cat_one]
                        for kw in kw_list:
                            if kw not in main_dict["keywords"][main_id]:
                                main_dict["keywords"][main_id].append(kw)
                FR.close()

        
        #--> glycation
        for sheet_name in ["proteoform_glycation_sites_"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "glycation"
                prop_dict = {
                    "start_pos":"glycation_site_uniprotkb",
                    "end_pos":"glycation_site_uniprotkb",
                    "glytoucan_ac":"saccharide",
                    "type":"glycation_type",
                    "residue":"amino_acid",
                    "comment":"uniprotkb_glycation_annotation_comment"

                }
                xref_info = ["xref_key", "xref_id"]
                combo_flist_one = ["glycation_site_uniprotkb","saccharide","glycation_type",
                    "amino_acid"
                ]
                combo_flist_two = ["xref_key", "xref_id"]
                load_obj = main_dict[prop_name]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)
                
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id,pos,glytoucan_ac,gly_type,amino_acid = combo_id.split("|")
                    if pos.strip() == "":
                        continue

                    cond_list = []
                    cond_list.append(pos == "")
                    cond_list.append(amino_acid == "")
                    cond_list.append(gly_type == "")
                    if True in cond_list:
                        load_obj.pop(combo_id)
                    load_obj[combo_id]["start_pos"]  = int(pos)
                    load_obj[combo_id]["end_pos"]  = int(pos)
                    gly_type = gly_type[0].upper() + gly_type[1:].lower()
                    amino_acid = amino_acid[0].upper() + amino_acid[1:]
                    load_obj[combo_id]["type"] = gly_type
                    load_obj[combo_id]["residue"] = amino_acid
                    load_obj[combo_id]["relation"] = "attached"

                    if main_id not in main_dict["keywords"]:
                        main_dict["keywords"][main_id] = []

                    seen_kw = {"glycatedprotein":True}
                    # Add glycatedprotein_reported if evidence in these resources
                    for o in load_obj[combo_id]["evidence"]:
                        if o["database"] in ["UniCarbKB", "PubMed", "PDB"]:
                            seen_kw["glycatedprotein_reported"] = True
                            break
                    # Add "glycatedprotein_reported" and 
                    # "glycatedprotein_reported_with_glycan"
                    # if it has glycan attached
                    if load_obj[combo_id]["glytoucan_ac"] != "":
                        seen_kw["glycatedprotein_reported"] = True
                        seen_kw["glycatedprotein_reported_with_glycan"] = True
                    if "glycatedprotein_reported" not in seen_kw:
                        seen_kw["glycatedprotein_predicted"] = True

                    for kw in seen_kw:
                        if kw not in main_dict["keywords"][main_id]:
                            main_dict["keywords"][main_id].append(kw)      


        #--> pathway
        for sheet_name in ["protein_xref_kegg", "protein_xref_reactome"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "pathway"
                load_obj = main_dict[prop_name]
                prop_dict = {"id":"xref_id", "name":"xref_label"}
                combo_flist_one = ["xref_key", "xref_id"]
                xref_info = []
                combo_flist_two = []
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)

                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    xref_key, xref_id = combo_id.split("|")[-2], combo_id.split("|")[-1]
                    xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    load_obj[combo_id]["url"] = xref_url
                    load_obj[combo_id]["resource"] = xref_badge
                    if xref_id in pathway_dict:
                        load_obj[combo_id]["description"] = pathway_dict[xref_id]["description"]
                        #consider reactions in which main_id participates as an enzyme
                        for rxn_id in pathway_dict[xref_id]["reaction_list"]:
                            if rxn_id in reaction_dict:
                                if "enzyme_list" in reaction_dict[rxn_id]:
                                    if main_id in reaction_dict[rxn_id]["enzyme_list"]:
                                        if "reaction_list" not in load_obj[combo_id]:
                                            load_obj[combo_id]["reaction_list"] = []
                                        load_obj[combo_id]["reaction_list"].append(rxn_id)



        #--> pro_annotation
        for sheet_name in ["protein_pro_annotation"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "pro_annotation"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        xref_key = tmp_row[tmp_fl.index("xref_key")]
                        xref_id = tmp_row[tmp_fl.index("xref_id")]
                        pro_name = tmp_row[tmp_fl.index("pro_protein_name")]
                        pro_def = tmp_row[tmp_fl.index("pro_protein_definition")]
                        combo_id = "%s|%s" % (main_id, xref_id)
                        xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                        xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                        ev_obj = {"id":xref_id,"database":xref_badge, "url":xref_url}
                        o = {"name":pro_name, "definition":pro_def, "evidence":[ev_obj]}
                        load_obj[combo_id] = o


        #--> go_annotation
        for sheet_name in ["protein_go_annotation"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    uniprotkb_ac = main_dict["uniprot_ac"][main_id]
                    go_ann_dict = {}
                    seen_go_id = {}
                    for tmp_row in data_frame["data"][main_id]:
                        go_term_id = tmp_row[tmp_fl.index("go_term_id")]
                        go_term_label = tmp_row[tmp_fl.index("go_term_label")]
                        go_term_category = tmp_row[tmp_fl.index("go_term_category")]
                        pmid = tmp_row[tmp_fl.index("pmid")]
                        url = libgly.get_xref_url(map_dict, "protein_xref_uniprotkb", uniprotkb_ac,is_cited)
                        ev_obj = {"database": "UniProtKB", "id":uniprotkb_ac, "url":url}
                        go_term_id = go_term_id.replace("_", ":")
                        url = libgly.get_xref_url(map_dict,"protein_xref_go", go_term_id,is_cited)
                        o = {"name": go_term_label, "id":go_term_id, "url":url, "evidence":[ev_obj],
                                "pmid":pmid}
                        if go_term_category not in go_ann_dict:
                            go_ann_dict[go_term_category] = []
                        if go_term_category not in seen_go_id:
                            seen_go_id[go_term_category] = {}

                        if go_term_id not in seen_go_id[go_term_category]:
                            go_ann_dict[go_term_category].append(o)
                            seen_go_id[go_term_category][go_term_id] = True

                    for go_term_category in go_ann_dict:
                        parts = go_term_category.split("_")
                        cat_label = parts[0][0:1].upper() + parts[0][1:]
                        cat_label += " " + parts[1][0:1].upper() + parts[1][1:]
                        n = len(go_ann_dict[go_term_category])
                        obj = {"name":cat_label, "total":n, "go_terms":[]}
                        for o in go_ann_dict[go_term_category]:
                            obj["go_terms"].append(o)
                        if main_id not in main_dict["go_annotation"]:
                            main_dict["go_annotation"][main_id] = {"categories":[]}
                        main_dict["go_annotation"][main_id]["categories"].append(obj)
                        combo_list = ["uniprotkb_canonical_ac", "go_term_label"]
                        update_record_stat(record_stat,  file_name, "go_annotation", 1, combo_list)

        #--> pathways
        for sheet_name in ["protein_pathways_reactome"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "pathways"
                load_obj = main_dict[prop_name]
                prop_dict = {"start_pos":"start_pos", "end_pos":"end_pos"}
                combo_flist_one = ["start_pos", "end_pos"]
                xref_info = []
                combo_flist_two = []
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)

                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    load_obj[combo_id]["annotation"] = "n_glycosylation_sequon"


        #--> site_annotation
        for sheet_name in ["protein_glycosylation_motifs"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "site_annotation"
                load_obj = main_dict[prop_name]
                prop_dict = {"start_pos":"start_pos", "end_pos":"end_pos"}
                combo_flist_one = ["start_pos", "end_pos"]
                xref_info = []
                combo_flist_two = []
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)

                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    load_obj[combo_id]["annotation"] = "n_glycosylation_sequon"
                    load_obj[combo_id]["start_pos"]  = int(load_obj[combo_id]["start_pos"])
                    load_obj[combo_id]["end_pos"]  = int(load_obj[combo_id]["end_pos"])

        #--> ptm_annotation
        for sheet_name in ["protein_ptm_annotation"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "ptm_annotation"
                load_obj = main_dict[prop_name]
                prop_dict = {"annotation":"ptm_annotation"}
                combo_flist_one = ["ptm_annotation"]
                xref_info = ["xref_key", "xref_id"]
                combo_flist_two = ["xref_key", "xref_id"]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)


        #--> mutagenesis
        for sheet_name in ["protein_site_annotation"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "mutagenesis"
                load_obj = main_dict[prop_name]
                prop_dict = {
                    "start_pos":"start_pos", "end_pos":"end_pos", 
                    "sequence_org":"ref_aa", "sequence_mut":"alt_aa",
                    "ann_type":"ann_type",
                    "comment":"annotation"
                }
                combo_flist_one = ["start_pos", "end_pos", "ref_aa", "alt_aa","annotation",
                        "ann_type"]
                xref_info = ["xref_key", "xref_id"]
                combo_flist_two = ["xref_key", "xref_id"]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)
        
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    if load_obj[combo_id]["ann_type"] not in ["Mutagenesis_Annotation"]:
                        load_obj[combo_id]["deleteflag"] = True
                    #load_obj[combo_id].pop("ann_type")
                    load_obj[combo_id]["start_pos"] = int(load_obj[combo_id]["start_pos"])
                    load_obj[combo_id]["end_pos"] = int(load_obj[combo_id]["end_pos"])
        
        
        #--> snv
        for sheet_name in ["protein_mutation_cancer", "protein_mutation_somatic",
                "protein_mutation_germline"]:
            if file_name.find(sheet_name) != -1 and file_name.find("glycoeffect") == -1:
                species = file_name.split("_")[0]
                prop_name = "snv"
                prop_dict_one = {
                    "start_pos":"aa_pos",
                    "end_pos":"aa_pos",
                    "subjects_tested":"patients_tested",
                    "subjects_positive":"patients_positive",
                    "frequency":"mut_freq",
                    "sequence_org":"ref_aa",
                    "sequence_mut":"alt_aa",
                    "comment":"filter_flags",
                    "chr_id":"chr_id",
                    "chr_pos":"chr_start_pos",
                    "ref_nt":"ref_nt",
                    "alt_nt":"alt_nt",
                    "minor_allelic_frequency":"minor_allelic_frequency"
                }
                prop_dict_two = {
                    "start_pos":"begin_aa_pos",
                    "end_pos":"end_aa_pos",
                    "sequence_org":"ref_aa",
                    "sequence_mut":"alt_aa",
                    "comment":"filter_flags",
                    "chr_id":"chr_id",
                    "chr_pos":"chr_pos",
                    "ref_nt":"ref_nt",
                    "alt_nt":"alt_nt",
                    "minor_allelic_frequency":"minor_allelic_frequency"
                } 
                prop_dict = prop_dict_one if sheet_name.find("cancer") != -1 else prop_dict_two
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        val_dict = {}
                        for p in prop_dict:
                            val_dict[p] = tmp_row[tmp_fl.index(prop_dict[p])]
                        #ignore rows with empty chr_id
                        if "chr_id" in prop_dict:
                            if val_dict["chr_id"].strip() == "":
                                continue
                        combo_list_one = [main_id]
                        combo_list_two = [main_id]
                        for p in ["start_pos", "sequence_org", "sequence_mut"]:
                            combo_list_one.append(val_dict[p])
                            combo_list_two.append(val_dict[p])
                        for k in ["xref_key", "xref_id"]:
                            combo_list_two.append(tmp_row[tmp_fl.index(k)])
                        combo_id_one = "|".join(combo_list_one)
                        combo_id_two = "|".join(combo_list_two)

                        if combo_id_one not in load_obj["seen"]:
                            load_obj["seen"][combo_id_one] = True
                            load_obj[combo_id_one] = {"evidence":[], "glycoeffect":[]}
                            for p in val_dict:
                                load_obj[combo_id_one][p] = val_dict[p]
                        
                        if combo_id_two not in load_obj["seen"]:
                            load_obj["seen"][combo_id_two] = True
                            xref_key = combo_list_two[-2]
                            xref_id = combo_list_two[-1]
                            if xref_key not in map_dict["xrefkey2url"]:
                                continue
                            xref_url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                            xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                            o = {"id":xref_id,"database":xref_badge, "url":xref_url}
                            load_obj[combo_id_one]["evidence"].append(o)
                        combo_id = combo_id_one
                        ann_score = len(load_obj[combo_id]["comment"].split(";"))
                        load_obj[combo_id]["ann_score"] = ann_score
                        load_obj[combo_id]["start_pos"] = int(load_obj[combo_id]["start_pos"])
                        load_obj[combo_id]["end_pos"] = int(load_obj[combo_id]["end_pos"])
                        load_obj[combo_id]["site_lbl"] = "%s%s" % (load_obj[combo_id]["sequence_org"], load_obj[combo_id]["start_pos"])

                        do_id = tmp_row[tmp_fl.index("do_id")]
                        mim_id = tmp_row[tmp_fl.index("mim_id")]
                        mondo_id =""
                        somatic_status = tmp_row[tmp_fl.index("somatic_status")]
                        
                        kw_list = ["ebi"]
                        if sheet_name.find("cancer") != -1:
                            kw_list = ["biomuta", "somatic", "cancer"]
                        elif somatic_status == "1" and "somatic" not in kw_list:
                            kw_list.append("somatic")
                        elif "germline" not in kw_list:
                            kw_list.append("germline")
                        disease_obj = get_disease_info(do_id, mondo_id, mim_id)
                        if disease_obj != {}:
                            kw_list.append("disease")
                            if "disease" not in load_obj[combo_id]:
                                load_obj[combo_id]["disease"] = [disease_obj]
                            else:
                                d_url_list = []
                                for o in load_obj[combo_id]["disease"]:
                                    d_url = o["recommended_name"]["url"]
                                    d_url_list.append(d_url)
                                if do_id != "":
                                    xref_key,xref_id = "protein_xref_do", do_id
                                    url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                                    if url not in d_url_list:
                                        load_obj[combo_id]["disease"].append(disease_obj)
                                elif mim_id != "":
                                    xref_key,xref_id = "protein_xref_omim", mim_id
                                    url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                                    if url not in d_url_list:
                                        load_obj[combo_id]["disease"].append(disease_obj)

                        if "keywords" not in load_obj[combo_id]:
                            load_obj[combo_id]["keywords"] = []
                        for kw in kw_list:
                            if kw not in load_obj[combo_id]["keywords"]:
                                load_obj[combo_id]["keywords"].append(kw)
                
                data_frame = {}
                in_file = in_file.replace(".csv", "_glycoeffect.csv")
                libgly.load_sheet_as_dict(data_frame,in_file,",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        aa_pos = tmp_row[tmp_fl.index("aa_pos")]
                        ref_aa = tmp_row[tmp_fl.index("ref_aa")]
                        alt_aa = tmp_row[tmp_fl.index("alt_aa")]
                        effect = tmp_row[tmp_fl.index("effect")]
                        combo_id = "%s|%s|%s|%s" % (main_id,aa_pos, ref_aa, alt_aa)
                        if combo_id not in load_obj:
                            continue 
                        if effect not in load_obj[combo_id]["glycoeffect"]:
                            load_obj[combo_id]["glycoeffect"].append(effect)




        #--> crossref
        for sheet_name in ["_protein_xref_"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "crossref"
                load_obj = main_dict[prop_name]
                prop_dict = {"id":"xref_id"}
                combo_flist_one = ["xref_key", "xref_id"]
                xref_info = []
                combo_flist_two = []
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)

                combo_id_list = list(load_obj.keys())
                for combo_id in combo_id_list:
                    if combo_id == "seen":
                        continue
                    main_id = combo_id.split("|")[0]
                    uniprotkb_ac = main_id.split("-")[0]
                    xref_key, xref_id = combo_id.split("|")[-2], combo_id.split("|")[-1]
                    xref_badge = libgly.get_xref_badge(map_dict, xref_key)
                    if xref_badge.strip() == "":
                        load_obj.pop(combo_id)
                        continue 
                    xref_url = ""
                    if xref_key in ["protein_xref_brenda"]:
                        xref_url = map_dict["xrefkey2url"][xref_key][0]
                        xref_url = xref_url % (xref_id, uniprotkb_ac)
                    else:
                        xref_url =  libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                    
                    load_obj[combo_id]["url"] = xref_url
                    load_obj[combo_id]["database"] = xref_badge
                    cats_dict = map_dict["xrefkey2category"]
                    xref_categories = cats_dict[xref_key] if xref_key in cats_dict else []
                    load_obj[combo_id]["categories"] = xref_categories

 

        #--> biomarkers
        for sheet_name in ["protein_biomarkers"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "biomarkers"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]: 
                    #print (main_id, main_id in biomarker_dict)
                    if main_id in biomarker_dict:
                        for doc in biomarker_dict[main_id]:
                            combo_id = "%s|%s" % (main_id, doc["biomarker_id"])
                            if "aclist" in doc:
                                doc.pop("aclist")
                            load_obj[combo_id] = doc



        #-->expression_tissue
        for sheet_name in ["protein_expression_normal"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "expression_tissue"
                load_obj = main_dict[prop_name]
                prop_dict = {
                    "tissue":"uberon_anatomical_name",
                    "score":"expression_score",
                    "present":"expression_level_gene_relative"
                }
                combo_flist_one = ["uberon_anatomical_id", "expression_level_gene_relative"]
                xref_info = ["xref_key", "xref_id"]
                combo_flist_two = ["xref_key", "xref_id"]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)

                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    tissue_id = combo_id.split("|")[-2].split(":")[1]
                    tissue_name = load_obj[combo_id]["tissue"]
                    if tissue_id in map_dict["uberonid2name"]:
                        tissue_name = map_dict["uberonid2name"][tissue_id][0]
                    tissue_name = tissue_name[0].upper() + tissue_name[1:]
                    name_space = "UBERON"
                    t_xref_key = "tissue_xref_" + name_space.lower()
                    t_xref_url = libgly.get_xref_url(map_dict, t_xref_key, tissue_id,is_cited)
                    o = {
                        "name":tissue_name,
                        "namespace":name_space,
                        "id":tissue_id, 
                        "url":t_xref_url
                    }
                    load_obj[combo_id]["tissue"] = o
                    #call = "yes" if load_obj[combo_id]["present"] == "present" else "no"
                    call = load_obj[combo_id]["present"] 
                    load_obj[combo_id]["present"] = call
                
                    if tissue_id in map_dict["uberonid2doid"]:
                        for do_id in map_dict["uberonid2doid"][tissue_id]:
                            load_obj["seen"][do_id] = True




        #--> expression_disease
        for sheet_name in ["protein_expression_disease"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "expression_disease"
                load_obj = main_dict[prop_name]
                prop_dict = {
                    "significant":"significance", 
                    "trend":"direction",
                    "do_id":"do_id",
                    "parent_doid":"parent_doid"
                }
                combo_flist_one = ["do_id"]
                xref_info = ["xref_key", "xref_id"]
                combo_flist_two = ["xref_key", "xref_id"]
                n = load_property_objlist(load_obj, in_file, prop_dict,xref_info,
                    combo_flist_one, combo_flist_two)
                combo_list = ["uniprotkb_canonical_ac"] + combo_flist_one
                update_record_stat(record_stat,  file_name, prop_name, n, combo_list)
                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    do_id = load_obj[combo_id]["do_id"]
                    parent_doid = load_obj[combo_id]["parent_doid"]
                    url = libgly.get_xref_url(map_dict, "protein_xref_do", do_id,is_cited)
                    mondo_id, mim_id = "", ""
                    disease_obj = get_disease_info(do_id, mondo_id, mim_id)
                    if disease_obj != {}:
                        if "disease" not in load_obj[combo_id]:
                            load_obj[combo_id]["disease"] = [disease_obj]
                        else:
                            d_url_list = []
                            for o in load_obj[combo_id]["disease"]:
                                d_url = o["recommended_name"]["url"]
                                d_url_list.append(d_url)
                            if do_id != "":
                                xref_key,xref_id = "protein_xref_do", do_id
                                url = libgly.get_xref_url(map_dict, xref_key, xref_id,is_cited)
                                if url not in d_url_list:
                                    load_obj[combo_id]["disease"].append(disease_obj)

                    tissue_seen_obj = main_dict["expression_tissue"]["seen"]
                    cond_list = [do_id in tissue_seen_obj]
                    cond_list += [parent_doid in tissue_seen_obj]
                    load_obj[combo_id].pop("do_id")
                    load_obj[combo_id].pop("parent_doid")
                    if cond_list == [False, False]:
                        load_obj[combo_id]["deleteflag"] = True




    tmp_dict = {}
    for main_id in main_dict["uniprot_canonical_ac"]:
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



    sec_list = main_dict.keys()
    for sec in sec_list:
        main_dict[sec].pop("seen")
        for combo_id in sorted(main_dict[sec]):
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
        if "refseq" in obj:
            if obj["refseq"] == {} and main_id in refseq_dict:
                obj["refseq"] = refseq_dict[main_id]




        obj["snv"].sort(key=get_sort_key_value_mut, reverse=True)
        for o in obj["snv"]:
            o.pop("ann_score")
        for o in obj["mutagenesis"]:
            o.pop("ann_type")

        obj["synthesized_glycans"] = []
        if main_id in synthesized_glycans:
            obj["synthesized_glycans"] = synthesized_glycans[main_id]
      
        obj["enzyme_annotation"] = enzyme_ann[main_id] if main_id in enzyme_ann else []
 
        obj["cluster_types"] = []
        seen_cluster_types = {}
        if "isoforms" in obj:
            if len(obj["isoforms"]) > 1:
                tmp_o = {"name":"isoformset.uniprotkb","label":"Isoforms"}
                if tmp_o["name"] not in seen_cluster_types:
                    obj["cluster_types"].append(tmp_o)
                    seen_cluster_types[tmp_o["name"]] = True

        if "orthologs" in obj:
            for o in obj["orthologs"]:
                for oo in o["evidence"]:
                    lbl = oo["database"] + " Homologs"
                    name = "homologset." + oo["database"].lower()
                    tmp_o = {"name":name,"label":lbl}
                    if tmp_o["name"] not in seen_cluster_types:
                        obj["cluster_types"].append(tmp_o)
                        seen_cluster_types[tmp_o["name"]] = True

        obj["sequence_features"] = get_seq_features(obj)
    
        section_stats.get_protein_sec_stats(obj, "protein")

        out_file = "jsondb/proteindb/%s.json" % (main_id)
        out_str = json.dumps(obj, indent=4)
        if len(out_str) > 16000000:
            out_file = "jsondb/jumbodb/proteindb/%s.json" % (main_id)

        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        record_count += 1 
    msg = "make-proteindb: final %s protein objects for %s" % (record_count, species_name)
    csvutil.write_log_msg(log_file, msg, "a")





    return


def load_pdbid2url(pdbid2url):

    file_list = glob.glob("downloads/pdb/current/*.pdb")
    for in_file in file_list:
        pdb_id = in_file.split("/")[-1][:-4]
        pdbid2url[pdb_id] = "https://data.glygen.org/ln2downloads/pdb/current/%s.pdb" % (pdb_id)
    return



#######################################
def main():


    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-s","--sec",action="store",dest="sec",help="Object section (OPTIONAL)")
    (options,args) = parser.parse_args()

                        
    global config_obj   
    global species_obj  
    global species_list
    global map_dict
    global doid2xrefid
    global xrefid2doid
    global data_dir
    global main_dict
    global record_stat
    global disease_id2names
    global disease_id2desc
    global recname_dict
    global submittedname_dict
    global disease_id2names
    global disease_id2desc
    global DEBUG
    global is_cited
    global pdbid2url
    global glycan_dict
    global glycan_class
    global synthesized_glycans
    global seq_dict
    global header_dict    
    global custom_xref_key_dict
    global tool_url_dict

    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    data_dir = wrk_dir

    custom_xref_key_dict = {
        "arabidopsis":"protein_xref_ensembl_arabidopsis",
        "fruitfly":"protein_xref_ensembl_fruitfly",
        "sarscov2":"protein_xref_ensembl_sarscov2",
        "yeast":"protein_xref_ensembl_yeast"
    }   
    tool_url_dict = {
        "isoglyp":"https://github.com/jonmohl/ISOGlyP",
        "GlycoSiteMiner":"https://github.com/glygener/glycositeminer"
    }
 
    DEBUG = False
    #DEBUG = True



    log_file = "logs/make-proteindb.log"
    msg = "make-proteindb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")


    msg = "make-proteindb: loading pdbid2url"
    csvutil.write_log_msg(log_file, msg, "a")
    pdbid2url = {}
    load_pdbid2url(pdbid2url)

    species_obj, species_list = {}, []
    load_species_info(species_obj, species_list)
   
    #species = "bovine"
    #for nm in ["common_name", "glygen_name"]:
    #    name_val = species_obj[str(species_obj[species]["tax_id"])][nm]
    #    print (nm, name_val)
    #exit()             
    

    msg = "make-proteindb: loading seqinfo"
    csvutil.write_log_msg(log_file, msg, "a")
    seq_dict, header_dict = {}, {}
    load_sequence_info(seq_dict, header_dict)

    glycan_dict = {}
    msg = "make-proteindb: loading glycan list"
    csvutil.write_log_msg(log_file, msg, "a")
    load_glycan_masterlist(glycan_dict)

    glycan_class = {}
    msg = "make-proteindb: loading glycan_class"
    csvutil.write_log_msg(log_file, msg, "a")
    load_glycan_class(glycan_class)

    synthesized_glycans = {}
    msg = "make-proteindb: loading synthesized_glycans"
    csvutil.write_log_msg(log_file, msg, "a")
    load_synthesized_glycans(synthesized_glycans, glycan_class)

    is_cited = libgly.get_is_cited() if DEBUG == False else {}
    if DEBUG:
        species_list = ["human"]

    #species_list = ["yeast"] 

    for species in species_list:
        map_dict = {}
        main_dict = {}
        record_stat = {}
        disease_id2names, disease_id2desc = {}, {}
    
        msg = "make-proteindb: now processing species=%s" % (species)
        csvutil.write_log_msg(log_file, msg, "a")

        make_objects(species)


    return


if __name__ == '__main__':
    main()



