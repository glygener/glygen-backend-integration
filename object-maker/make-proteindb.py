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



##################
def get_sort_key_value_pub(obj):
    return obj["date"]

def get_sort_key_value_mut(obj):
    return obj["ann_score"]

def load_pathways(pathway_dict):

    file_list = glob.glob("reviewed/*_protein_pathways_reactome.csv")
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



def load_enzymes(enzyme_dict):

    canon2gene_name = {}
    file_list = glob.glob("reviewed/*_protein_masterlist.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            canon2gene_name[canon] = row[f_list.index("gene_name")]

    file_list = glob.glob("reviewed/*_protein_enzyme_annotation_uniprotkb.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            e_id = canon.split("-")[0]
            xref_key, xref_id = "protein_xref_uniprotkb", e_id
            xref_url =  map_dict["xrefkey2url"][xref_key][0] % (xref_id)
            xref_badge = map_dict["xrefkey2badge"][xref_key][0]
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
            
    return


def load_reactions(reaction_dict, participant_dict, ac2rxn, enzyme_dict):
    
    role_dict = {}
    file_list = glob.glob("reviewed/*_protein_participants_reactome.csv")
    file_list += glob.glob("reviewed/*_protein_participants_rhea.csv")

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




    file_list = glob.glob("reviewed/*_protein_reactions_reactome.csv")
    file_list += glob.glob("reviewed/*_protein_reactions_rhea.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        source = in_file.split("_")[-1].split(".")[0]
        for row in data_frame["data"]:
            rxn_id = row[f_list.index("reaction_id")]
            xref_id = rxn_id
            xref_key = "protein_xref_%s" % (source)

            xref_url =  map_dict["xrefkey2url"][xref_key][0] % (xref_id)
            xref_badge = map_dict["xrefkey2badge"][xref_key][0]
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

def load_protein_names():

    file_list = glob.glob("reviewed/*_protein_*names.csv")
    file_list += glob.glob("reviewed/*_protein_*names_refseq.csv")
    file_list += glob.glob("reviewed/*_protein_*names_uniprotkb.csv")
    sheet_info = load_name_sheet_info()
    for in_file in file_list:
        sheet_name = "protein_" + in_file.split("_protein_")[1].replace(".csv", "")
        prop_name = "gene_names" if sheet_name.find("genenames") != -1 else "protein_names"
        load_obj = main_dict[prop_name]
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
                xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                name_obj = { "name":name, "resource":xref_badge,  "type":name_type,
                        "url":xref_url, "id":xref_id}
                combo_id = "%s|%s|%s|%s" % (canon,name_type,name,resource)
                load_obj[combo_id] = name_obj


    return



def get_sorting_key(obj):
    return obj['sortorder']



def load_species_info(species_obj, species_list):

    seen = {}
    in_file = path_obj["misc"]+ "/species_info.csv"
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


def load_dictionaries(map_dict, misc_dir):

    dict_list_obj = json.loads(open("conf/protein_dictionaries.json", "r").read())
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


def load_disease_names():

    data_frame = {}
    in_file = data_dir + "/protein_disease_names.csv"
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
    in_file = data_dir + "/protein_disease_idmap.csv"
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
    in_file = data_dir + "/human_protein_genomics_england_disease.csv"
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
                if in_file.find("_glycosylation_sites_") != -1 and f == "saccharide":
                    if val not in seen_glycan:
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
                        obj_one["url"] = map_dict["xrefkey2url"]["protein_xref_do_placeholder"][0]
                    else:
                        do_name = ""
                        if do_id in disease_id2names["do"]:
                            do_name = disease_id2names["do"][do_id]["name"]
                            do_name = do_name[0].upper() + do_name[1:]+" [DO disease name]"
                        obj_one["name"] = do_name
                        obj_one["url"] = map_dict["xrefkey2url"]["protein_xref_do"][0] % (do_id)

            if combo_id not in seen_dict:
                seen_dict[combo_id] = True
                tmp_obj_dict[combo_id] = obj_one
                if combo_flist_two != []:
                    obj_one["evidence"] = []


            if combo_flist_two != []:
                xref_key = tmp_row[tmp_fl.index(xref_info[0])]
                xref_id = tmp_row[tmp_fl.index(xref_info[1])]
                if xref_key in map_dict["xrefkey2badge"]:
                    xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                    xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
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
  
    s = "*"
    file_list = glob.glob(data_dir + "%s_protein_allsequences.fasta" % (s))
    for in_file in file_list:
        for record in SeqIO.parse(in_file, "fasta"):
            seq_id = record.id.split("|")[1]
            seq_dict[seq_id] = str(record.seq.upper())
    
    file_list = glob.glob(data_dir + "%s_protein_sequenceinfo.csv" % (s))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            isoform = row[f_list.index("uniprotkb_isoform_ac")]
            header_dict[isoform] = row[f_list.index("sequence_header")]

    return



def load_isoformlocus_info(isoformlocus_dict):
    
    s = species_list[0] if len(species_list) == 1 else "*"
    file_list = glob.glob(data_dir + "%s_protein_transcriptlocus.csv" % (s))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        tmp_fl = data_frame["fields"]
        for tmp_row in data_frame["data"]:
            isoform = tmp_row[tmp_fl.index("uniprotkb_isoform_ac")]
            transcript_id = tmp_row[tmp_fl.index("transcript_id")]
            peptide_id = tmp_row[tmp_fl.index("peptide_id")]
            t_url = map_dict["xrefkey2url"]["protein_xref_ensembl"][0] % (transcript_id)
            p_url = map_dict["xrefkey2url"]["protein_xref_ensembl"][0] % (peptide_id)
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

def load_glycan_masterlist():

    is_motif = {}
    in_file = path_obj["reviewed"] +  "glycan_motif.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("motif_ac_xref")]
        is_motif[ac] = True


    glycan_list = []
    data_frame = {}
    in_file = path_obj["reviewed"] +  "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        #if ac not in glycan_list and ac not in is_motif:
        if ac not in glycan_list:
            glycan_list.append(ac)

    return glycan_list


def load_genelocus_info(genelocus_dict):

    s = species_list[0] if len(species_list) == 1 else "*"
    file_list = glob.glob(data_dir + "%s_protein_genelocus.csv" % (s))                

    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        tmp_fl = data_frame["fields"]
        for tmp_row in data_frame["data"]:
            canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
            ensemble_gene_id = tmp_row[tmp_fl.index("ensembl_gene_id")]
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
            gene_url = map_dict["xrefkey2url"]["protein_xref_ensembl"][0] % (ensemble_gene_id)
            if canon not in genelocus_dict:
                genelocus_dict[canon] = {
                    "chromosome":chr_id
                    ,"start_pos":start_pos
                    ,"end_pos":end_pos
                    ,"strand":strand_sign
                    ,"evidence":[
                        {
                            "database":"Ensembl Gene"
                            ,"id":ensemble_gene_id
                            ,"url":gene_url
                        }
                    ]
                }
    return




def load_cluster_info(cls_dict, seq_dict):

    canon2genename = {}
    for in_file in glob.glob(data_dir + "/*_protein_masterlist.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        tmp_fl = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[tmp_fl.index("uniprotkb_canonical_ac")]
            canon2genename[canon] = row[tmp_fl.index("gene_name")]

    homolog_dict = {}
    member_dict = {}
    data_frame = {}
    in_file = data_dir + "/protein_homolog_clusters.csv"
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
        prop = "protein_names"
        recname = ""
        if canon in recname_dict["protein"]:
            recname =  recname_dict["protein"][canon]
        elif canon in submittedname_dict["protein"]:
            recname = submittedname_dict["protein"][canon]
        tax_id = row[tmp_fl.index("tax_id")]
        tax_name = species_obj[str(tax_id)]["long_name"]
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
            common_name = species_obj[str(tax_id)]["common_name"]
            if common_name != "":
                homolog_dict[canon]["common_name"] = common_name
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
                xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                o = {"id":xref_id, "url":xref_url, "database":xref_badge}
                main_dict["orthologs"][combo_id]["evidence"].append(o)


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
        name_xref_url = map_dict["xrefkey2url"][name_xref_key][0] % (name_xref_id)
        name_xref_badge = map_dict["xrefkey2badge"][name_xref_key][0]
        name = disease_id2names["do"][name_xref_id]["name"] if name_xref_id in disease_id2names["do"] else ""
        disease_desc = disease_id2desc["do"][name_xref_id] if name_xref_id in disease_id2desc["do"] else ""
        disease_info = {
                "key":name_xref_key,"id":"DOID:"+name_xref_id, "badge":name_xref_badge, 
            "url":name_xref_url, "name":name, "description":disease_desc
        }
    elif mondo_id != "":
        name_xref_key = "protein_xref_mondo"
        name_xref_id = mondo_id
        name_xref_url = map_dict["xrefkey2url"][name_xref_key][0] % (name_xref_id)
        name_xref_badge = map_dict["xrefkey2badge"][name_xref_key][0]
        name = disease_id2names["mondo"][name_xref_id]["name"]
        disease_desc = disease_id2desc["mondo"][name_xref_id]
        disease_info = {
                "key":name_xref_key,"id":"MONDO:"+name_xref_id, "badge":name_xref_badge,
            "url":name_xref_url, "name":name, "description":disease_desc
        }
    elif mim_id != "" and mim_id in disease_id2names["omim"]:
        name_xref_key = "protein_xref_omim"
        name_xref_id = mim_id
        name_xref_url = map_dict["xrefkey2url"][name_xref_key][0] % (name_xref_id)
        name_xref_badge = map_dict["xrefkey2badge"][name_xref_key][0]
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
            syn_xref_url = map_dict["xrefkey2url"][syn_xref_key][0] % (syn_xref_id)
            syn_xref_badge = map_dict["xrefkey2badge"][syn_xref_key][0]
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
    global species_list
    global map_dict
    global doid2xrefid
    global xrefid2doid
    global data_dir
    global misc_dir
    global main_dict
    global seen_glycan
    global disease_id2names
    global disease_id2desc
    global recname_dict
    global submittedname_dict


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"
    misc_dir = "generated/misc/"



    DEBUG = True if sec_name_list != [] else False


    glycan_list = load_glycan_masterlist()


    file_name_list = []
    ds_obj_list = json.loads(open(misc_dir + "/dataset-masterlist.json", "r").read())
    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if ds_name in ["homolog_alignments", "isoform_alignments"]:
            continue
        if obj["categories"]["species"] == []:
            if obj["integration_status"]["status"] == "integrate_all":
                if "protein" in obj["target_objects"]:
                    file_name_list.append("%s_%s" % (mol, ds_name))
        elif obj["integration_status"]["status"] != "integrate_none":
            sp_list_one = sorted(obj["categories"]["species"])
            for species in sp_list_one:
                if species in ["fruitfly"]:
                    continue
                if species not in obj["integration_status"]["excludelist"]:
                    if "protein" in obj["target_objects"]:
                        file_name_list.append("%s_%s_%s" % (species, mol, ds_name))

    #Since *_protein_expression_disease depends on *_protein_expression_normal,
    #move it to last order in the file_name_list
    tmp_list = []
    for f in file_name_list:
        if f.find("protein_expression_disease") != -1:
            tmp_list.append(f)
            file_name_list.remove(f)
    file_name_list += tmp_list



    sec_info = json.loads(open("generated/misc/protein_sectioninfo.json", "r").read())

    species_obj, species_list = {}, []
    load_species_info(species_obj, species_list)

    pattern_list = []
    if sec_name_list != []:
        for sec_name in sec_name_list:
            pattern_list += sec_info[sec_name]["sheetlist"]
        species_list = ["human"]
    else:
        for sec in sec_info:
            pattern_list += sec_info[sec]["sheetlist"]
    pattern_list = list(set(pattern_list))


    selected_file_name_list = []
    for file_name in file_name_list:
        cond_list = []
        for pat in ["_protein_masterlist"] + pattern_list:
            cond_list += [file_name.find(pat) != -1]
        if list(set(cond_list)) != [False]:
            selected_file_name_list.append(file_name)

    #Check missing files
    mising_files = []
    for file_name in selected_file_name_list:
        file_ext = "fasta" if file_name.find("protein_allsequences") != -1 else "csv"
        in_file = "%s%s.%s" % (data_dir,file_name,file_ext)
        if os.path.isfile(in_file) == False:
        #if os.path.isfile(in_file) == False and in_file.find("fruitfly") == -1:
            mising_files.append(in_file)

    
    if mising_files != []:
        print ("The following files are missing:")
        print (json.dumps(mising_files, indent=4))
        exit()

    ###############
    # Only for debug
    #tmp_list = []
    #for f in selected_file_name_list:
    #    if f.find("sarscov2_proteoform_glycosylation_sites_unicarbkb") != -1:
    #        tmp_list.append(f)
    #selected_file_name_list = tmp_list
    #print (json.dumps(selected_file_name_list, indent=4))
    #exit()

    ###############


    map_dict = {}
    load_dictionaries(map_dict, misc_dir)



    main_dict = {}
    record_stat = {}
    for sec in sec_info:
        main_dict[sec] = {"seen":{}}

    enzyme_dict = {}
    load_enzymes(enzyme_dict)
    kw = "enzyme"
    for main_id in enzyme_dict:
        if main_id not in main_dict["keywords"]:
            main_dict["keywords"][main_id] = []
            if kw not in main_dict["keywords"][main_id]:
                main_dict["keywords"][main_id].append(kw)
                combo_list = ["uniprotkb_canonical_ac", "keywords"]
                update_record_stat(record_stat,  file_name, "keywords", 1, combo_list)


    reaction_dict, participant_dict, ac2rxn = {}, {}, {}
    load_reactions(reaction_dict,participant_dict, ac2rxn, enzyme_dict)



    pathway_dict = {}
    load_pathways(pathway_dict)



    disease_id2names = {}
    disease_id2desc = {}
    load_disease_names()

    
    seq_dict, header_dict = {}, {}
    load_sequence_info(seq_dict, header_dict)

    genelocus_dict = {}
    load_genelocus_info(genelocus_dict)

    isoformlocus_dict = {}
    load_isoformlocus_info(isoformlocus_dict)

    recname_dict = {"protein":{}, "gene":{}}
    submittedname_dict = {"protein":{}}
    load_protein_names()


    refseq_dict = {}
    for combo_id in main_dict["protein_names"]:
        canon = combo_id.split("|")[0]
        if combo_id.split("|")[-1] == "refseq":
            obj = main_dict["protein_names"][combo_id]
            if canon not in refseq_dict:
                ac = obj["id"]
                url = map_dict["xrefkey2url"]["protein_xref_refseq"][0] % (ac)
                refseq_dict[canon] = {"ac":ac, "name":obj["name"], "url":url}



    


    cls_dict = {}
    load_cluster_info(cls_dict, seq_dict)


    seen_glycan = {}
    in_file = data_dir + "/glycan_masterlist.csv"
    sheet_obj = {}
    libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "glytoucan_ac")
    tmp_fl = sheet_obj["fields"]
    for glytoucan_ac in sheet_obj["data"]:
        if glytoucan_ac not in seen_glycan:
            seen_glycan[glytoucan_ac] = True

    data_frame = {}
    in_file = data_dir + "/glycan_classification.csv"
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
    tmp_fl = data_frame["fields"]
    glycan_class = {}
    for main_id in data_frame["data"]:
        for tmp_row in data_frame["data"][main_id]:
            g_type = tmp_row[tmp_fl.index("glycan_type")].strip()
            g_subtype = tmp_row[tmp_fl.index("glycan_subtype")].strip()
            g_type = "Other" if g_type == "" else g_type
            g_subtype = "Other" if g_subtype == "" else g_subtype
            glycan_class[main_id] = {"type":g_type, "subtype":g_subtype}


    synthesized_glycans = {}
    seen_glytoucan = {}
    in_file = data_dir + "/glycan_enzyme.csv"
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



    #--> is_reported dict for glycosylation
    is_reported = {}
    predicted_xref_key_list = [ "protein_xref_uniprotkb_gly"]
    for file_name in selected_file_name_list:
        file_ext = "fasta" if file_name.find("protein_allsequences") != -1 else "csv"
        in_file = "%s%s.%s" % (data_dir,file_name,file_ext)
        if file_name.find("_glycosylation_sites_") != -1:
            data_frame = {}
            libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
            tmp_fl = data_frame["fields"]
            for main_id in data_frame["data"]:
                for tmp_row in data_frame["data"][main_id]:
                    aa_pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                    glytoucan_ac = tmp_row[tmp_fl.index("saccharide")]
                    start_pos = tmp_row[tmp_fl.index("start_pos")]
                    end_pos = tmp_row[tmp_fl.index("end_pos")]
                    xref_key = tmp_row[tmp_fl.index("xref_key")]
                    if xref_key not in predicted_xref_key_list:
                        combo_id = "%s|%s|%s|%s" % (main_id,start_pos, end_pos,glytoucan_ac)
                        is_reported[combo_id] = True
                        combo_id = "%s|%s|%s|%s" % (main_id,start_pos, end_pos,"")
                        is_reported[combo_id] = True



    seen_combo_glycosylation = {}
    seen_combo_disease = {}

    file_idx = 1
    file_count = len(selected_file_name_list)
    
    for file_name in selected_file_name_list:
        
        file_ext = "fasta" if file_name.find("protein_allsequences") != -1 else "csv" 
        in_file = "%s%s.%s" % (data_dir,file_name,file_ext)
        if os.path.isfile(in_file) == False:
            print ("make-proteindb: file %s does NOT exist!" % (in_file))
            sys.exit()
  
        print ("make-proteindb: %s [%s/%s]" % (in_file, file_idx, file_count))
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
                gene_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
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


                url = map_dict["xrefkey2url"]["protein_xref_uniprotkb"][0] 
                url = url % (uniprotkb_ac)
                o = {
                    "name":species_obj[species]["long_name"]
                    ,"taxid":species_obj[species]["tax_id"]
                    ,"evidence":[
                        {"database":"UniProtKB", "id":uniprotkb_ac,"url":url}
                    ]   
                }
                common_name = species_obj[str(species_obj[species]["tax_id"])]["common_name"]
                if common_name != "":
                    o["common_name"] = common_name
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
                    refseq_url = map_dict["xrefkey2url"]["protein_xref_refseq"][0] % (refseq_ac)
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
                    xref_url =  map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                    load_obj[combo_id]["date"] = load_obj[combo_id]["date"].split(" ")[0]
                    xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                    load_obj[combo_id]["reference"] = [
                        {
                            "type":xref_badge,
                            "id":xref_id,
                            "url":xref_url
                        }
                    ]
      
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
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
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
                #for combo_id in load_obj:
                #    if combo_id == "seen":
                #        continue
                #    xref_key,xref_id = combo_id.split("|")[-2], combo_id.split("|")[-1]
                #    ann_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                #    load_obj[combo_id]["url"] = ann_url

        #--> phosphorylation
        for sheet_name in ["_phosphorylation_sites_"]:
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
                        load_obj[combo_id]["start_pos"]  = int(load_obj[combo_id]["start_pos"])
                        load_obj[combo_id]["end_pos"]  = int(load_obj[combo_id]["end_pos"])




         
        for sheet_name in ["_glycosylation_sites_"]:
            if file_name.find(sheet_name) != -1:
                species = file_name.split("_")[0]
                prop_name = "glycosylation"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        #src_file_name = tmp_row[tmp_fl.index("src_file_name")]
                        tmp_aa_pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                        tmp_glytoucan_ac = tmp_row[tmp_fl.index("saccharide")]
                        tmp_amino_acid = tmp_row[tmp_fl.index("amino_acid")]
                        site_seq = tmp_row[tmp_fl.index("site_seq")]
                        start_pos = tmp_row[tmp_fl.index("start_pos")]
                        end_pos = tmp_row[tmp_fl.index("end_pos")]
                        start_aa = tmp_row[tmp_fl.index("start_aa")]
                        end_aa = tmp_row[tmp_fl.index("end_aa")]
                        aa_pos_list = tmp_aa_pos.split("|")
                        glytoucan_ac_list = tmp_glytoucan_ac.split("|")
                        amino_acid_list = tmp_amino_acid.split("|")
                        aa_pos = aa_pos_list[0]
                        glytoucan_ac = glytoucan_ac_list[0]
                        amino_acid = amino_acid_list[0]
                        
                        if glytoucan_ac != "" and glytoucan_ac not in glycan_list:
                            continue

                        gly_type = tmp_row[tmp_fl.index("glycosylation_type")]
                        gly_type = gly_type[0].upper() + gly_type[1:]
                        gly_subtype = ""
                        if "glycosylation_subtype" in tmp_fl:
                            gly_subtype = tmp_row[tmp_fl.index("glycosylation_subtype")]
                        
                        xref_key = tmp_row[tmp_fl.index("xref_key")]
                        xref_id = tmp_row[tmp_fl.index("xref_id")]
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                        parts = tmp_row[tmp_fl.index("src_xref_key")].split("_")
                        src_suffix = "_".join(parts[2:])
                        parts = tmp_row[tmp_fl.index("xref_key")].split("_")
                        xref_suffix = "_".join(parts[2:])
                        combo_id = "%s|%s|%s|%s" % (main_id,start_pos, end_pos,glytoucan_ac)
                        

                        site_cat_one, site_cat_two = "", ""
                        if combo_id in is_reported:
                            #site_cat_one = "reported.%s" % (src_suffix)
                            site_cat_one = "reported"
                            site_cat_two = "reported"
                            if src_suffix == "automatic_literature_mining":
                                site_cat_two = "automatic_literature_mining"
                            
                            if glytoucan_ac != "":
                                site_cat_one = "reported_with_glycan"
                                #site_cat_one += ".%s" % (src_suffix)
                                site_cat_two += "_with_glycan"
                        else:
                            #site_cat_one = "predicted.%s" % (src_suffix)
                            site_cat_one = "predicted"
                            site_cat_two = "predicted"
                            if glytoucan_ac != "":
                                site_cat_one = "predicted_with_glycan"
                                #site_cat_one += ".%s" % (src_suffix)
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
                   
                        if combo_id_one not in seen_combo_glycosylation:
                            gly_obj = {
                                "glytoucan_ac":glytoucan_ac,
                                "type":gly_type,
                                "subtype":gly_subtype,
                                "site_category":site_cat_two,
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
                        ev_obj = {"id":xref_id, "database":xref_badge, "url":xref_url}
                        combo_id_two = "%s|%s" % (xref_key,xref_id)
                        
                        if combo_id_two not in seen_combo_glycosylation[combo_id_one]:
                            load_obj[combo_id_one]["evidence"].append(ev_obj)
                            seen_combo_glycosylation[combo_id_one].append(combo_id_two)

                        if main_id not in main_dict["keywords"]:
                            main_dict["keywords"][main_id] = []
                        kw_list = ["glycoprotein", "glycoprotein_" + site_cat_one]
                        for kw in kw_list:
                            if kw not in main_dict["keywords"][main_id]:
                                main_dict["keywords"][main_id].append(kw)


        
        #--> glycation
        for sheet_name in ["_glycation_sites_"]:
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
                    xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                    xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                    load_obj[combo_id]["url"] = xref_url
                    load_obj[combo_id]["resource"] = xref_badge
                    if xref_id in pathway_dict:
                        load_obj[combo_id]["description"] = pathway_dict[xref_id]["description"]
                        #consider reactions in which main_id participates as an enzyme
                        for rxn_id in pathway_dict[xref_id]["reaction_list"]:
                            if rxn_id in reaction_dict:
                                if "enzyme_list" in reaction_dict[rxn_id]:
                                    if main_id in reaction_dict[rxn_id]["enzyme_list"]:
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
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
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
                        url = map_dict["xrefkey2url"]["protein_xref_uniprotkb"][0] % (uniprotkb_ac)
                        ev_obj = {"database": "UniProtKB", "id":uniprotkb_ac, "url":url}
                        go_term_id = go_term_id.replace("_", ":")
                        url = map_dict["xrefkey2url"]["protein_xref_go"][0] % (go_term_id)
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
                    "chr_pos":"chr_pos",
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
                            xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                            xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                            o = {"id":xref_id,"database":xref_badge, "url":xref_url}
                            load_obj[combo_id_one]["evidence"].append(o)

                        combo_id = combo_id_one
                        ann_score = len(load_obj[combo_id]["comment"].split(";"))
                        load_obj[combo_id]["ann_score"] = ann_score
                        load_obj[combo_id]["start_pos"] = int(load_obj[combo_id]["start_pos"])
                        load_obj[combo_id]["end_pos"] = int(load_obj[combo_id]["end_pos"])

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
                                    url = map_dict["xrefkey2url"][xref_key][0]%(xref_id)
                                    if url not in d_url_list:
                                        load_obj[combo_id]["disease"].append(disease_obj)
                                elif mim_id != "":
                                    xref_key,xref_id = "protein_xref_omim", mim_id
                                    url = map_dict["xrefkey2url"][xref_key][0]%(xref_id)
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

                for combo_id in load_obj:
                    if combo_id == "seen":
                        continue
                    main_id = combo_id.split("|")[0]
                    uniprotkb_ac = main_id.split("-")[0]
                    xref_key, xref_id = combo_id.split("|")[-2], combo_id.split("|")[-1]
                    xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                    xref_url = map_dict["xrefkey2url"][xref_key][0]
                    if xref_url.find("%s") != -1:
                        if xref_key in ["protein_xref_brenda"]:
                            xref_url = xref_url % (xref_id, uniprotkb_ac)
                        else:
                            xref_url = xref_url % (xref_id)
                    load_obj[combo_id]["url"] = xref_url
                    load_obj[combo_id]["database"] = xref_badge
                

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
                    uberon_id = combo_id.split("|")[-2].split(":")[1]
                    uberon_name = load_obj[combo_id]["tissue"]
                    if uberon_id in map_dict["uberonid2name"]:
                        uberon_name = map_dict["uberonid2name"][uberon_id][0]
                    uberon_name = uberon_name[0].upper() + uberon_name[1:]
                    url = map_dict["xrefkey2url"]["protein_xref_uberon"][0] % (uberon_id)
                    o = {"name":uberon_name,"uberon":uberon_id, "url":url}
                    load_obj[combo_id]["tissue"] = o
                    #call = "yes" if load_obj[combo_id]["present"] == "present" else "no"
                    call = load_obj[combo_id]["present"] 
                    load_obj[combo_id]["present"] = call
                
                    if uberon_id in map_dict["uberonid2doid"]:
                        for do_id in map_dict["uberonid2doid"][uberon_id]:
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
                    url = map_dict["xrefkey2url"]["protein_xref_do"][0] % (do_id)
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
                                url = map_dict["xrefkey2url"][xref_key][0]%(xref_id)
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



    sec_list = main_dict.keys() if DEBUG == False else sec_name_list
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


        out_file = path_obj["jsondbpath"] + "/proteindb/%s.json" % (main_id)
        obj["snv"].sort(key=get_sort_key_value_mut, reverse=True)
        for o in obj["snv"]:
            o.pop("ann_score")
        for o in obj["mutagenesis"]:
            o.pop("ann_type")

        obj["synthesized_glycans"] = []
        if main_id in synthesized_glycans:
            obj["synthesized_glycans"] = synthesized_glycans[main_id]
       
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


        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        record_count += 1 
    print ("make-proteindb: final filtered in: %s protein objects" % (record_count))


    out_file = path_obj["jsondbpath"] + "/logs/proteindb.json" 
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(record_stat, indent=4)))





if __name__ == '__main__':
        main()



