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
        if mol == "proteoform":
            continue
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
    

    
    
    sec_info = json.loads(open("generated/misc/motif_sectioninfo.json", "r").read())
    
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


    motif_glytoucan_ac_list = load_motif_glytoucan_ac_list()


    record_stat = {}
    glycan_type_dict = {}
    

    motif_name_dict = {}

    
    main_obj_dict = {}
    for file_name in file_name_list:
        
        cond_list = []
        for pat in ["_protein_masterlist"] + pattern_list:
            cond_list += [file_name.find(pat) != -1]
        if DEBUG and list(set(cond_list)) == [False]:
            continue

        
        in_file = "%s%s.csv" % (data_dir,file_name)
        if os.path.isfile(in_file) == False:
            print "make-motifdb: file %s does NOT exist!" % (in_file)
            sys.exit()
  
        print "make-motifdb:", in_file

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
                if number_monosaccharides != "":
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
                            url = map_dict["xrefkey2url"]["glycan_xref_inchi_key"][0]
                            url = url % (inchi_key)
                            o = {"key":inchi_key, "url":url}
                            main_dict["inchi_key"][main_id] = o
                            combo_list = ["glytoucan_ac", field]
                            update_record_stat(record_stat, file_name, "inchi_key", 1, combo_list)
        
        

        #--> publication
        #--> keep this empty for now
        for sheet_name in ["glycan_citations_glycomotif"]:
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
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (main_id)
                        main_dict["glycans"][combo_id] = {"glytoucan_ac":main_id, 
                                "url":xref_url
                        }
                        if motif_ac_xref not in main_dict["motif_ac"]:
                            main_dict["motif_ac"][motif_ac_xref] = []
                        if motif_ac not in main_dict["motif_ac"][motif_ac_xref]:
                            main_dict["motif_ac"][motif_ac_xref].append(motif_ac)
                        xref_url = map_dict["xrefkey2url"][xref_key][0] % (motif_ac_xref)
                        xref_badge = map_dict["xrefkey2badge"][xref_key][0]
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
                                
                                            xref_url = map_dict["xrefkey2url"][xref_key][0] % (xref_id)
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
        #--> don't do crossref here
        for sheet_name in ["glycan_xrefffff_"]:
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
        if main_id not in motif_glytoucan_ac_list:
            continue
        obj = tmp_dict[main_id]
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

            out_file = path_obj["jsondbpath"] + "/motifdb/%s.json" % (motif_ac)
            for k in ["mass", "mass_pme", "number_monosaccharides"]:
                if obj[k] == 0.0 or obj[k] == 0:
                    obj.pop(k)
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(obj, indent=4)))
            record_count += 1 
    print "make-motifdb: final filtered in: %s motif objects" % (record_count)
    
    out_file = path_obj["jsondbpath"] + "/logs/motifdb.json" 
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(record_stat, indent=4)))





if __name__ == '__main__':
        main()



