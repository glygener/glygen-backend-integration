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

    

def load_motif_glytoucan_ac_list():
    
    
    seen = {}
    data_frame = {}
    in_file = path_obj["reviewed"]+ "/glycan_motif.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        motif_ac_xref = row[f_list.index("motif_ac_xref")]
        seen[motif_ac_xref] = True


    return seen.keys()


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

    dict_list_obj = json.loads(open("../../conf/glycan_dictionaries.json", "r").read())
    for dict_name in dict_list_obj:
        map_dict[dict_name] = {}
        ind_list = dict_list_obj[dict_name]["indexlist"]
        for pattern in dict_list_obj[dict_name]["fileglob"]:
            for in_file in glob.glob(misc_dir + pattern):
                print "Now processing ... %s" % (in_file)
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
                        obj_one["url"] = map_dict["xrefkey2url"]["protein_xref_do_placeholder"][0]
                    else:
                        do_name = doid2name[do_id][0]
                        obj_one["name"] = do_name[0].upper() + do_name[1:] + " [DO disease name]"
                        obj_one["url"] = map_dict["xrefkey2url"]["protein_xref_do"][0] % (do_id)
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
                xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                xref_url = map_dict["xrefkey2url"][xref_key][0]
                if map_dict["xrefkey2url"][xref_key][0].find("%s") != -1:
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



def load_protein_objects(file_name_list):

    seen = {}
    for file_name in file_name_list:
        for patt in ["proteoform_glycosylation_sites", "protein_matrixdb", "glycan_enzyme"]:
            if file_name.find(patt) != -1:
                in_file = data_dir + "/%s.csv" % (file_name)
                data_frame = {}
                print "Now processing ... %s" % (in_file)
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    seen[main_id] = True
    
    canon_list = seen.keys()

    protein_obj_dict = {}
    for canon in canon_list:
        protein_jsonfile = "jsondb/proteindb/%s.json" % (canon)
        if os.path.isfile(protein_jsonfile) == False:
            protein_obj_dict[canon] = {}
        else:
            protein_obj_dict[canon] = json.loads(open(protein_jsonfile,"r").read())


    return protein_obj_dict


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


    DEBUG = True
    #DEBUG = False

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
                if "glycan" in obj["target_objects"]:
                    file_name_list.append("%s_%s" % (mol, ds_name))
        else:
            sp_list_one = sorted(obj["categories"]["species"])
            for species in sp_list_one:
                if species not in obj["integration_status"]["excludelist"]:
                    if "glycan" in obj["target_objects"]:
                        file_name_list.append("%s_%s_%s" % (species, mol, ds_name))
    
    
    sec_info = json.loads(open("generated/misc/glycan_sectioninfo.json", "r").read())
    main_dict = {}
    for sec in sec_info:
        main_dict[sec] = {"seen":{}}



    species_obj, species_list = {}, []
    load_species_info(species_obj, species_list)

    map_dict = {}
    load_dictionaries(map_dict, misc_dir)



    pattern_list = []
    if sec_name != None:
        pattern_list += sec_info[sec_name]["sheetlist"]
    else:
        for sec in sec_info:
            pattern_list += sec_info[sec]["sheetlist"]
    pattern_list = list(set(pattern_list))


    fully_determined_list = []
    data_frame = {}
    in_file = "reviewed/glycan_fully_determined.csv" 
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
    tmp_fl = data_frame["fields"]
    for main_id in data_frame["data"]:
        fully_determined_list.append(main_id)
        

    selected_file_name_list = []
    for file_name in file_name_list:
        cond_list = []
        for pat in ["_protein_masterlist"] + pattern_list:
            cond_list += [file_name.find(pat) != -1]
        if list(set(cond_list)) != [False]:
            selected_file_name_list.append(file_name)

    #print json.dumps(selected_file_name_list, indent=4)
    #exit()

    in_file = "%s%s" % (data_dir,"glycan_subsumption.csv")
    heirarchy_dict = load_subsumption_heirarchy(in_file)



    protein_obj_dict = load_protein_objects(selected_file_name_list)
    motif_glytoucan_ac_list = load_motif_glytoucan_ac_list()

    record_stat = {}
    
    
    main_obj_dict = {}
    for file_name in selected_file_name_list:
        
        in_file = "%s%s.csv" % (data_dir,file_name)
        if os.path.isfile(in_file) == False:
            print "file %s does NOT exist!" % (in_file)
            sys.exit()
  
        print "Now processing ", in_file

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
                number_monosaccharides = tmp_row[tmp_fl.index("monosaccharides")]
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
                if number_monosaccharides != "":
                    prop_name = "number_monosaccharides"
                    number_monosaccharides = int(tmp_row[tmp_fl.index("monosaccharides")])
                    main_dict[prop_name][main_id] = number_monosaccharides
                    combo_list = ["glytoucan_ac", "number_monosaccharides"]
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
                "byonic":"sequence_byonic"
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
                            url = map_dict["xrefkey2url"]["glycan_xref_inchi_key"][0]
                            url = url % (inchi_key)
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
                    #load_obj[combo_id].pop("annotation_category")
                    if load_obj[combo_id]["taxid"]  == "":
                        load_obj[combo_id]["taxid"] = -1
                    else:
                        common_name = species_obj[load_obj[combo_id]["taxid"]]["common_name"]
                        if common_name != "":
                            load_obj[combo_id]["common_name"] = common_name
                        load_obj[combo_id]["taxid"] = int(load_obj[combo_id]["taxid"])


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
                        rec_name = get_protein_info(protein_obj_dict[canon], "recommendedname")
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
                    xref_url =  map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                    load_obj[combo_id]["date"] = load_obj[combo_id]["date"].split(" ")[0]
                    xref_badge = map_dict["xrefkey2badge"][xref_key]
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
                prop_dict = {"id":"motif_ac", "name":"motif_name", "synonym":"alternative_name"}
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
                xref_info = ["xref_key", "xref_id"]
                combo_flist_two = ["xref_key", "xref_id"]
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
                        #load_obj[combo_id]["deleteflag"] = True
                        #continue
                    protein_obj = protein_obj_dict[canon]
                    if protein_obj == {}:
                        load_obj[combo_id]["deleteflag"] = True
                        continue
                    gene_name = get_protein_info(protein_obj, "gene_name")
                    rec_name = get_protein_info(protein_obj, "recommendedname")
                    tax_id = get_protein_info(protein_obj, "tax_id")
                    tax_name = get_protein_info(protein_obj, "tax_name")
                    load_obj[combo_id]["protein_name"] = rec_name
                    load_obj[combo_id]["gene_name"] = gene_name
                    load_obj[combo_id]["tax_id"] = tax_id
                    load_obj[combo_id]["tax_name"] = tax_name
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
                load_obj = main_dict[prop_name]
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
                        cond_one = main_id != related_accession
                        cond_two = glytoucan_type == relationship
                        cond_three = glytoucan_type in gt_list
                        if cond_one and cond_two and cond_three:
                            combo_id = "%s|%s|%s" % (main_id, related_accession,relationship)
                            rl = ""
                            if main_id in heirarchy_dict:
                                if related_accession in heirarchy_dict[main_id]:
                                    rl = heirarchy_dict[main_id][related_accession]
                            load_obj[combo_id] = {
                                "id":related_accession, 
                                "type":glytoucan_type, 
                                "relationship":rl
                            }


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
                prop_dict = {"uniprot_canonical_ac":"uniprotkb_canonical_ac"}
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
                    protein_obj = protein_obj_dict[canon]
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
                    xref_key = tmp_row[tmp_fl.index("xref_key")]
                    xref_id = tmp_row[tmp_fl.index("xref_id")]
                    xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                    xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                    ev_obj = {"id":xref_id, "database":xref_badge, "url":xref_url}
                    load_obj[combo_id]["evidence"] = [ev_obj]

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
                    uniprotkb_ac = main_id.split("-")[0]
                    xref_key, xref_id = combo_id.split("|")[-2], combo_id.split("|")[-1]
                    xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                    xref_url = map_dict["xrefkey2url"][xref_key][0]
                    if xref_url.find("%s") != -1:
                        xref_url = xref_url % (xref_id)
                    load_obj[combo_id]["url"] = xref_url
                    load_obj[combo_id]["database"] = xref_badge
        

        #--> expression
        seen_evdn = {}
        for sheet_name in ["proteoform_glycosylation_sites_glyconnect", 
                "proteoform_glycosylation_sites_unicarbkb"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "expression"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "saccharide")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    for tmp_row in data_frame["data"][main_id]:
                        if main_id == "":
                            continue
                        canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                        
                        tmp_aa_pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                        tmp_aa_three = tmp_row[tmp_fl.index("amino_acid")]
                        aa_pos_list = tmp_aa_pos.split("|")
                        aa_three_list = tmp_aa_three.split("|")
                        aa_pos = aa_pos_list[0]
                        aa_three = aa_three_list[0]
                        if aa_pos == "":
                            continue
                        abundance, tissue_name, uberon_id = "", "", ""
                        if "source_tissue_uberon_id" in tmp_fl:
                            tissue_name = tmp_row[tmp_fl.index("source_tissue_name")]
                            uberon_id = tmp_row[tmp_fl.index("source_tissue_uberon_id")]
                        
                        if "abundance" in tmp_fl:
                            abundance = tmp_row[tmp_fl.index("abundance")]
                        cl_name, cl_id = "", ""
                        if "source_cell_line_cellosaurus_id" in tmp_fl:
                            cl_name = tmp_row[tmp_fl.index("source_cell_line_name")]
                            cl_id = tmp_row[tmp_fl.index("source_cell_line_cellosaurus_id")]
                        if "cellosaurus_id" in tmp_fl:
                            cl_name = tmp_row[tmp_fl.index("cellosuaurus_cell_line_name")]
                            cl_id = tmp_row[tmp_fl.index("cellosaurus_id")]
                        src_xref_key = tmp_row[tmp_fl.index("src_xref_key")]
                        src_xref_id = tmp_row[tmp_fl.index("src_xref_id")]
                        if src_xref_key == "protein_xref_glyconnect":
                            src_xref_key = "glycan_xref_glyconnect"
                            src_xref_id = tmp_row[tmp_fl.index("structure_id")]
                        src_xref_badge = map_dict["xrefkey2badge"][src_xref_key][0]
                        src_xref_url = map_dict["xrefkey2url"][src_xref_key][0] % (src_xref_id)
                        xref_key = tmp_row[tmp_fl.index("xref_key")]
                        xref_id = tmp_row[tmp_fl.index("xref_id")]
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
                        combo_id = "%s|%s|%s|%s" % (main_id,canon,aa_pos,uberon_id)
                        if combo_id not in load_obj["seen"]:
                            cl_key = "glycan_xref_cellosaurus"
                            cl_url = map_dict["xrefkey2url"][cl_key][0] % (cl_id)
                            cl_obj = {"name": cl_name,"cellosaurus_id":cl_id,"url":cl_url}
                            cl_obj = cl_obj if cl_id not in ["", "0"] else {}
                            uberon_key = "glycan_xref_uberon"
                            uberon_url = map_dict["xrefkey2url"][uberon_key][0] % (uberon_id)
                            tissue_obj = { "name":tissue_name,"uberon_id":uberon_id,
                                "url":uberon_url
                            }
                            tissue_obj = tissue_obj if uberon_id != "" else {}

                            load_obj[combo_id] = {
                                "uniprot_canonical_ac":canon, 
                                "start_pos":int(aa_pos),
                                "end_pos":int(aa_pos),
                                "residue":aa_three,
                                "tissue":tissue_obj,
                                "cell_line":cl_obj,
                                "abundance":abundance,
                                "evidence":[
                                    {"id":src_xref_id, "database":src_xref_badge,
                                        "url":src_xref_url}
                                ]
                            }
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
                            evdn_combo = "%s|%s|%s|%s|%s|%s"%(main_id,canon,aa_pos,uberon_id,src_xref_key,src_xref_id)
                            load_obj["seen"][evdn_combo] = True

                        evdn_combo = "%s|%s|%s|%s|%s|%s"%(main_id,canon,aa_pos,uberon_id,xref_key,xref_id)
                        if evdn_combo not in load_obj["seen"]:
                            ev_obj = {"id":xref_id, "database":xref_badge, "url":xref_url}
                            load_obj[combo_id]["evidence"].append(ev_obj)
                            load_obj["seen"][evdn_combo] = True



        #--> composition               
        for sheet_name in ["glycan_monosaccharide_composition"]:
            if file_name.find(sheet_name) != -1:
                prop_name = "composition"
                load_obj = main_dict[prop_name]
                data_frame = {}
                libgly.load_sheet_as_dict(data_frame, in_file, ",", "glytoucan_ac")
                tmp_fl = data_frame["fields"]
                for main_id in data_frame["data"]:
                    type_combo_list = []
                    for tmp_row in data_frame["data"][main_id]:
                        for residue in map_dict["residue2name"]:
                            n = int(tmp_row[tmp_fl.index(residue)])
                            name = map_dict["residue2name"][residue][0]
                            cid = map_dict["residue2cid"][residue][0]
                            residue = "other" if residue.lower() == "xxx" else residue.lower()
                            url = map_dict["xrefkey2url"]["glycan_xref_monosaccharide_residue_name"][0] % (cid)
                            o = {"name":name, "residue":residue, "count":n}
                            if cid != "":
                                o = {"name":name, "residue":residue, "count":n, "cid":cid, "url":url}
                            combo_id = "%s|%s|%s" % (main_id, residue, cid)
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
                        t_tag = map_dict["subtypetags"][glycan_type][0]
                        type_url = map_dict["xrefkey2url"]["glycan_xref_glycan_type"][0] % (t_tag)
                        o = {"type":{"name":glycan_type, "url":type_url}}
                        if t_tag == "xxx":
                            o["type"].pop("url")
                        if glycan_subtype != "no subtype":
                            s_tag = map_dict["subtypetags"][glycan_subtype][0]
                            subtype_url = map_dict["xrefkey2url"]["glycan_xref_glycan_type"][0] % (s_tag)
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
                            t_tag = map_dict["subtypetags"][glycan_type][0]
                            type_url = map_dict["xrefkey2url"]["glycan_xref_glycan_type"][0] % (t_tag)
                            s_tag = map_dict["subtypetags"][glycan_subtype][0]
                            subtype_url = map_dict["xrefkey2url"]["glycan_xref_glycan_type"][0] % (s_tag)
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

        out_file = path_obj["jsondbpath"] + "/glycandb/%s.json" % (main_id)
        obj = tmp_dict[main_id]
        if obj["classification"] == []:
            o = {"type":{"name":"Other", "url":""}, "subtype":{"name":"Other", "url":""}}
            obj["classification"] = [o]
    
        for k in ["mass", "mass_pme", "number_monosaccharides"]:
            if obj[k] == 0.0 or obj[k] == 0:
                obj.pop(k)
        
        obj["tool_support"] = {"gnome":"no","sandbox":"no"}
        if "mass" in obj:
            obj["tool_support"]["gnome"] = "yes"
        for xobj in obj["crossref"]:
            if xobj["database"] == "SandBox":
                obj["tool_support"]["sandbox"] = "yes"
                break
        #print json.dumps(obj, indent=4)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        record_count += 1 
        #if record_count%1000 == 0:
        #    print " ... created %s objects" % (record_count)
    print " ... final filtered in: %s glycan objects" % (record_count)


    
    out_file = path_obj["jsondbpath"] + "/logs/glycandb.json" 
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(record_stat, indent=4)))





if __name__ == '__main__':
        main()



