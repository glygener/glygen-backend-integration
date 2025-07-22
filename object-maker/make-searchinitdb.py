import os,sys
import string
from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json

import libgly
import csvutil

__version__="1.0"
__status__ = "Dev"


def get_doc_list(coll, k_list):


    file_list = glob.glob("jsondb/" + coll.replace("c_", "") + "db/" + PATTERN)
    doc_list = []
    for f in file_list:
        doc = json.loads(open(f, "r").read())
        if coll == "c_network":
            new_doc = {}
            for obj in doc["recordlist"]:
                for k in k_list:
                    new_doc[k] = obj[k]
            doc_list.append(new_doc)
        elif k_list != []:
            new_doc = {}
            for k in k_list:
                if k in doc:
                    new_doc[k] = doc[k]
            doc_list.append(new_doc) 
        else:
            doc_list.append(doc)
        
    
    return doc_list




def get_field_limits(limits_info):

    k_list_dict = {
        "c_glycan":["glytoucan_ac","crossref","glycoprotein","motifs","enzyme","publication"],
        "c_protein":[
            "uniprot_ac","uniprot_canonical_ac","refseq","protein_names","gene_names",
            "go_annotation","glycosylation","sequence","pathway","publication"
        ]
    }
    
    res_obj = {}
    for cat in limits_info:
        res_obj[cat] = {}
        coll = ""
        coll = "c_protein" if cat == "proteinsearch" else coll
        coll = "c_glycan" if cat == "glycansearch" else coll
        if coll == "":
            for lbl in limits_info[cat]:
                if lbl not in res_obj[cat]:
                    res_obj[cat][lbl] = 1000
            continue
        doc_list = get_doc_list(coll, k_list_dict[coll])
        for doc in doc_list:
            for lbl in limits_info[cat]:
                if lbl not in res_obj[cat]:
                    res_obj[cat][lbl] = 0
                if limits_info[cat][lbl] == {}:
                    res_obj[cat][lbl] = 1000
                else:
                    for prop_lineage in limits_info[cat][lbl]:
                        prop_list = prop_lineage.split(".")
                        sec_type = limits_info[cat][lbl][prop_lineage]
                        val = ""
                        if sec_type == 1:
                            val = doc[prop_lineage]
                            if len(val) > res_obj[cat][lbl]:
                                res_obj[cat][lbl] = len(val)
                        elif sec_type == 2:
                            if  prop_list[1] in doc[prop_list[0]]:
                                val = doc[prop_list[0]][prop_list[1]]
                                if len(val) > res_obj[cat][lbl]:
                                    res_obj[cat][lbl] = len(val)
                        elif sec_type == 3:
                            for o in doc[prop_list[0]]:
                                val = o[prop_list[1]]
                                if len(val) > res_obj[cat][lbl]:
                                    res_obj[cat][lbl] = len(val)
                        elif sec_type == 4:
                            if prop_list[1] in doc[prop_list[0]]:
                                for obj in doc[prop_list[0]][prop_list[1]]:
                                    for o in obj[prop_list[2]]:
                                        val =  o[prop_list[3]]
                                        if len(val) > res_obj[cat][lbl]:
                                            res_obj[cat][lbl] = len(val)
    

    return res_obj


def get_disease_search_init(collection):

    res_obj = {}
  
    seen = {"organism":{}, "biomarker_types":{}}
    file_list = glob.glob("jsondb/" + collection.replace("c_", "") + "db/" + PATTERN)
    for f in file_list:
        doc = json.loads(open(f, "r").read())
        for sec in ["proteins", "glycans"]:
            for obj in doc[sec]:
                for o in obj["species"]:
                    if o["glygen_name"] != "":
                        org_obj_str = json.dumps({ "name": o["glygen_name"], "id":o["taxid"]})
                        seen["organism"][org_obj_str] = True
        if doc["biomarkers"] != []:
            for obj in doc["biomarkers"]:
                for o in obj["best_biomarker_role"]:
                    seen["biomarker_types"][o["role"]] = True


    res_obj["organism"] = []
    for org_obj_str in seen["organism"]:
        res_obj["organism"].append(json.loads(org_obj_str))

    res_obj["biomarker_types"] =  list(seen["biomarker_types"].keys())

    res_obj["simple_search_category"] = [
         {"id":"any", "display":"Any"}
        ,{"id":"organism", "display":"Organism"}
        ,{"id":"protein", "display":"Protein"}
        ,{"id":"glycan", "display":"Glycan"}
        ,{"id":"biomarker", "display":"Biomarker"}
    ]


    return res_obj


def get_site_search_init(collection):
    
    res_obj = {}
    seen = {"aa_type_list":{}}
    doc_list = get_doc_list(collection,["site_seq"])
    for obj in doc_list:
        if "site_seq" in obj:
            k = obj["site_seq"]
            if len(k) > 1 or len(k) == 0:
                continue
            label = aa_dict[k]["name"] if k in aa_dict else ""
            seen["aa_type_list"][k] = {"id":k, "label":label}

    for k_one in seen:
        res_obj[k_one] = []
        for k_two in seen[k_one]:
            res_obj[k_one].append(seen[k_one][k_two])


    res_obj["annotation_type_list"] = [
        {"id":"snv_flag", "label":"SNV"},
        {"id":"glycosylation_flag", "label":"Glycosylation"},
        {"id":"mutagenesis_flag","label":"Mutagenesis"},
        {"id":"phosphorylation_flag","label":"Phosphorylation"},
        {"id":"glycation_flag", "label":"Glycation"}
    ]
    return res_obj



def get_glycan_search_init(collection):

    species_map = json.loads(open("generated/misc/species_map.json", "r").read())
    tax_id_map = {}
    for tax_id in species_map:
        tax_id_map[tax_id] = species_map[tax_id]["ref_tax_id"]
  
 
    tax_id_list = []
    for species in species_list:
        tax_id = species_obj[species]["tax_id"]
        tax_id_list.append(tax_id)
        if tax_id in tax_id_map:
            tax_id_list.append(tax_id_map[tax_id])
    tax_id_list = sorted(list(set(tax_id_list)))

    res2cid = {}
    lines = open("generated/misc/monosaccharide_residue_name.csv", "r").read().split("\n")[1:]
    for line in lines:
        if line.strip() == "":
            continue
        line = line.replace("\"", "")
        row = line.split(",")
        res, cid  = row[0].strip(), row[2].strip()
        res2cid[res.lower()] = cid


    res_obj = {
        "glycan_mass":{
            "native":{ "name":"Native", "min":0, "max":0 }
            ,"permethylated":{ "name":"Permethylated", "min":0, "max":0 }
        }
        ,"number_monosaccharides":{}, 
        "organism":[], 
        "glycan_type":[],
        "id_namespace":[]
    }
    
    
    min_mass = 10000000.0
    max_mass = -min_mass
    min_mass_pme = 10000000.0
    max_mass_pme = -min_mass_pme

    min_monosaccharides = 100000
    max_monosaccharides = -min_monosaccharides

    seen = {"biomarker_types":{}, "glycan_type":{}, "organism":{}, "id_namespace":{}}
    seen_glygen_name = {}

    k_list = [
        "glytoucan_ac", "mass", "mass_pme", "number_monosaccharides","species","composition",
        "classification", "crossref", "biomarkers"
    ]
    comp_dict = {}
    doc_list = get_doc_list(collection, k_list)
    for obj in doc_list:
        glytoucan_ac = obj["glytoucan_ac"]
        if "mass" in obj:
            obj["mass"] = float(obj["mass"])
            if obj["mass"] > 0.0:
                min_mass = round(obj["mass"], 2) if obj["mass"] < min_mass else min_mass
                max_mass = round(obj["mass"], 2) if obj["mass"] > max_mass else max_mass
       
        if "mass_pme" in obj:
            obj["mass_pme"] = float(obj["mass_pme"])
            if obj["mass_pme"] > 0.0:
                min_mass_pme = round(obj["mass_pme"], 2) if obj["mass_pme"] < min_mass_pme else min_mass_pme
                max_mass_pme = round(obj["mass_pme"], 2) if obj["mass_pme"] > max_mass_pme else max_mass_pme


        if "biomarkers" in obj:
            for o in obj["biomarkers"]:
                for role in o["best_biomarker_role"]:
                    seen["biomarker_types"][role] = True
            

        if "number_monosaccharides" in obj:
            min_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] < min_monosaccharides else min_monosaccharides
            max_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] > max_monosaccharides else max_monosaccharides

        if "species" in obj:
            for o in obj["species"]:
                #if o["taxid"] not in tax_id_list:
                #    continue
                if o["glygen_name"] != "":
                    if o["glygen_name"] not in seen_glygen_name:
                        org_obj_str = json.dumps({ "name": o["glygen_name"], "id":o["taxid"]})
                        seen["organism"][org_obj_str] = True
                        seen_glygen_name[o["glygen_name"]] = True


        if "composition" in obj:
            for o in obj["composition"]:
                res, res_name, res_count = o["residue"], o["name"], o["count"] 
                cid = res2cid[res.lower()] if res.lower() in res2cid else "" 
                res_url = "https://pubchem.ncbi.nlm.nih.gov/compound/%s" % (cid)
                if res not in comp_dict:
                    comp_dict[res] = {"residue":res, "min":10000, "max":0, "url":res_url }
                if res_count < comp_dict[res]["min"]:
                    comp_dict[res]["min"] = res_count
                if res_count > comp_dict[res]["max"]:
                    comp_dict[res]["max"] = res_count

        if "classification" in obj:
            for o in obj["classification"]:
                type_name = o["type"]["name"]
                if type_name not in seen["glycan_type"]:
                    seen["glycan_type"][type_name] = {}
                if "subtype" in o:
                    subtype_name = o["subtype"]["name"]
                    subtype_name = "other" if subtype_name.lower() == "other" else subtype_name
                    if subtype_name not in seen["glycan_type"][type_name]:
                        seen["glycan_type"][type_name][subtype_name] = True 
        
        if "crossref" in obj:
            for o in obj["crossref"]:
                xref_badge, xref_id = o["database"], o["id"]
                seen["id_namespace"][xref_badge] = True


    res_obj["glycan_mass"]["native"]["min"] = min_mass
    res_obj["glycan_mass"]["native"]["max"] = max_mass
    res_obj["glycan_mass"]["permethylated"]["min"] = min_mass_pme
    res_obj["glycan_mass"]["permethylated"]["max"] = max_mass_pme
     
    res_obj["number_monosaccharides"]["min"] = min_monosaccharides
    res_obj["number_monosaccharides"]["max"] = max_monosaccharides 

    res_obj["id_namespace"] = list(seen["id_namespace"].keys())
    res_obj["biomarker_types"] =  list(seen["biomarker_types"].keys())

    for type_name in seen["glycan_type"]:
        o = {"name":type_name,"subtype":[]}
        for subtype_name in seen["glycan_type"][type_name]:
            o["subtype"].append(subtype_name)
        res_obj["glycan_type"].append(o)
    
    for org_obj_str in seen["organism"]:
        res_obj["organism"].append(json.loads(org_obj_str))


    res_obj["simple_search_category"] = [
         {"id":"any", "display":"Any"}
        ,{"id":"glycan", "display":"Glycan"}
        ,{"id":"protein", "display":"Protein"}
        ,{"id":"enzyme", "display":"Enzyme"}
        ,{"id":"organism", "display":"Organism"}
    ]
    
    res_obj["composition"] = []
    for res in comp_dict:
        res_obj["composition"].append(comp_dict[res])


    return res_obj


def get_idmapping_search_init():

    res_obj = {
        "protein":{"id":"protein", "label":"Protein", "namespace":[]}, 
        "glycan":{"id":"glycan", "label":"Glycan","namespace":[]}
    }
    record_type_list = list(res_obj.keys())

    for record_type in record_type_list:
        seen = {}
        coll = "c_" + record_type
        xref_count = {}
        doc_list = get_doc_list(coll, ["crossref"])
        for doc in doc_list:
            tmp_dict = {}
            for obj in doc["crossref"]:
                xref_badge, xref_id = obj["database"], obj["id"]
                if xref_id.strip() == "":
                    continue
                tmp_dict[xref_badge] = True

                if xref_badge not in xref_count:
                    xref_count[xref_badge] = {}
                if xref_id not in xref_count[xref_badge]:
                    xref_count[xref_badge][xref_id] = 0
                xref_count[xref_badge][xref_id] += 1

            d_list = list(tmp_dict.keys())

            for d_one in d_list:
                if d_one not in seen:
                    seen[d_one] = {}
                for d_two in d_list:
                    if d_one != d_two:
                        seen[d_one][d_two] = True
        


        id_list_dict = {}
        for xref_badge in xref_count:
            id_list = sorted(xref_count[xref_badge], key=xref_count[xref_badge].get, reverse=True)
            if len(id_list) > 5:
                id_list = id_list[:5]
            id_list_dict[xref_badge] = id_list

        res_obj[record_type]["namespace"] = []
        for d_one in seen:
            id_list = id_list_dict[d_one] if d_one in id_list_dict else []
            o = {"source":d_one, "targetlist":list(seen[d_one].keys()), 
                    "example_id_list":id_list}
            res_obj[record_type]["namespace"].append(o)
    

    return res_obj

def get_usecases_search_init():
    res_obj = {
        "species_to_glycohydrolases":{"organism":{}}
        ,"species_to_glycosyltransferases":{"organism":{}}
        ,"species_to_glycoproteins":{"organism":{}}
    }
    query_dict = {
    "species_to_glycosyltransferases":"glycosyltransferase-activity",
    "species_to_glycohydrolases":"glycohydrolase-activity",
    "species_to_glycoproteins":"glycoprotein"
    }

    evdn_lbl_dict = {
        "all_sites": "All sites",
        "all_reported_sites_with_without_glycans": "All reported sites (with or without Glycans)",
        "sites_reported_with_glycans": "Sites reported with Glycans",
        "sites_reported_without_glycans": "Sites reported without Glycans",
        "predicted_sites": "Predicted sites",
        "sites_detected_by_literature_mining": "Sites detected by literature mining"
    }
    


 
    #doc_list = get_doc_list("c_protein", ["uniprot_canonical_ac", "species", "keywords"])
    file_list = glob.glob("jsondb/proteindb/*.json")
    for in_file in file_list:
        protein_doc = json.load(open(in_file))
        list_file = "jsondb/listdb/" + in_file.split("/")[-1]
        if os.path.isfile(list_file) == False:
            continue
        doc = json.load(open(list_file))
        obj = {"tax_id":doc["tax_id"], "tax_name":doc["organism"]}
        obj["keywords"] = protein_doc["keywords"]
        if len(protein_doc["glycosylation"]) > 0:
            obj["keywords"].append("all_sites")  
        #waiting for clarification for all_reported_sites_with_without_glycans 
        if doc["reported_n_glycosites_with_glycan"] + doc["reported_o_glycosites_with_glycan"] > 0:
            obj["keywords"].append("sites_reported_with_glycans")
        elif doc["reported_n_glycosites"] + doc["reported_o_glycosites"] > 0:
            obj["keywords"].append("sites_reported_without_glycans")
        if doc["mined_glycosites"] > 0:
            obj["keywords"].append("sites_detected_by_literature_mining")
        if doc["predicted_glycosites"] > 0:
            obj["keywords"].append("predicted_sites")

        if "sites_reported_with_glycans" in obj["keywords"] or "sites_reported_without_glycans" in obj["keywords"]:
            obj["keywords"].append("all_reported_sites_with_without_glycans")

        for svc in query_dict:
            if query_dict[svc] not in obj["keywords"]:
                continue
            tax_id, tax_name = obj["tax_id"], obj["tax_name"]        
            if tax_id not in res_obj[svc]["organism"]:
                res_obj[svc]["organism"][tax_id] = {"id":tax_id, "name":tax_name}
            if svc == "species_to_glycoproteins":
                if "evidence_type" not in res_obj[svc]["organism"][tax_id]:
                    res_obj[svc]["organism"][tax_id]["evidence_type"] = {}
                for kw in obj["keywords"]:
                    if kw in evdn_lbl_dict:
                        res_obj[svc]["organism"][tax_id]["evidence_type"][kw] = True

    svc = "species_to_glycoproteins"
    for tax_id in res_obj[svc]["organism"]:
        o_list = []
        for kw in res_obj[svc]["organism"][tax_id]["evidence_type"]:
            o_list.append({"id":kw, "display":evdn_lbl_dict[kw]})
        res_obj[svc]["organism"][tax_id]["evidence_type"] = o_list
 

    return res_obj


def get_super_search_init():


    doc_list = get_doc_list("c_network", ["record_type", "record_id"])
    seen_one = {}
    for doc in doc_list:
        src_record_type, src_record_id = doc["record_type"], doc["record_id"]
        if src_record_type not in seen_one:
            seen_one[src_record_type] = {}
        seen_one[src_record_type][src_record_id] = True
    res_obj = {}
    for src_record_type in seen_one:
        res_obj[src_record_type] = len(list(seen_one[src_record_type].keys()))


    return res_obj

def get_biomarker_search_init(collection):


    k_list = ["biomarker_id", "best_biomarker_role", "biomarker_component"]
    doc_list = get_doc_list(collection, k_list)
    

    res_obj = {}
    seen = {}
    f_list = ["best_biomarker_role", "assessed_entity_type"]
    for f in f_list:
        seen[f] = {}

    for doc in doc_list:
        f = "best_biomarker_role"
        if f in doc:
            for obj in doc[f]:
                val = obj["role"]
                seen[f][val] = True
        
        
        f = "assessed_entity_type"
        if "biomarker_component" in doc:
            for obj in doc["biomarker_component"]:
                val = obj[f]
                seen[f][val] = True
       

    for f in f_list:
        res_obj[f] =  list(seen[f].keys())

    res_obj["simple_search_category"] = [
         {"id":"any", "display":"Any"}
        ,{"id":"biomarker", "display":"Biomarker"}
        ,{"id":"condition", "display":"Condition"}
    ]



    return res_obj

def get_protein_search_init(collection):

    tax_id_list = []
    for species in species_list: 
        tax_id_list.append(species_obj[species]["tax_id"]) 


    res_obj = {"protein_mass":{}, "organism":[], 
            "glycosylation_types":{}, "id_namespace":[]}
    min_mass = 10000000.0
    max_mass = -min_mass

    seen = {"protein_mass":{}, "organism":{}, "amino_acid":{}, 
            "glycosylation_types":{}, 
            "biomarker_types":{},
            "id_namespace":{}
    }
    k_list = ["uniprot_canonical_ac", "mass","species","glycosylation", "crossref", "biomarkers"]

    for tax_id in tax_id_list:
        q_obj = {"species.taxid":{"$eq":tax_id}}
        doc_list = get_doc_list(collection, k_list)
        for obj in doc_list:
            uniprot_canonical_ac = obj["uniprot_canonical_ac"]
            if "mass" in obj:
                mass = float(obj["mass"]["chemical_mass"])
                if mass == -1.0:
                    continue
                min_mass = round(mass, 2) if mass < min_mass else min_mass
                max_mass = round(mass, 2) if mass > max_mass else max_mass
            if "species" in obj:
                for o in obj["species"]:
                    if o["taxid"] in tax_id_list:
                        if o["glygen_name"] != "":
                            org_obj_str = json.dumps({ "name": o["glygen_name"], "id":o["taxid"]})
                            seen["organism"][org_obj_str] = True
            if "biomarkers" in obj:
                for o in obj["biomarkers"]:
                    for role in o["best_biomarker_role"]:
                        seen["biomarker_types"][role] = True

            if "glycosylation" in obj:
                for o in obj["glycosylation"]:
                    if "residue" in o:
                        seen["amino_acid"][o["residue"]] = True
                    if "type" in o:
                        if o["type"].find(";") == -1:
                            if o["type"] != "":
                                g_type = o["type"][0].upper() + o["type"][1:].lower()
                                if g_type not in seen["glycosylation_types"]:
                                    seen["glycosylation_types"][g_type] = {}
                                if "subtype" in o:
                                    if o["subtype"] != "":
                                        g_subtype = o["subtype"]
                                        g_subtype = "other" if g_subtype.lower() == "other" else g_subtype
                                        if g_subtype[0:2] == "O-" or g_subtype[0:2] == "C-":
                                            g_subtype = g_subtype[0:2] + g_subtype[2].upper() + g_subtype[3:]
                                        seen["glycosylation_types"][g_type][g_subtype] = True

            if "crossref" in obj:
                for o in obj["crossref"]:
                    xref_badge, xref_id = o["database"], o["id"]
                    seen["id_namespace"][xref_badge] = True


    res_obj["id_namespace"] = list(seen["id_namespace"].keys())
    res_obj["protein_mass"]["min"] = min_mass
    res_obj["protein_mass"]["max"] = max_mass
    for org_obj_str in seen["organism"]:
        res_obj["organism"].append(json.loads(org_obj_str))
    for g_type in seen["glycosylation_types"]:
        res_obj["glycosylation_types"][g_type] = list(seen["glycosylation_types"][g_type].keys())

    res_obj["biomarker_types"] =  list(seen["biomarker_types"].keys())

    res_obj["simple_search_category"] = [
         {"id":"any", "display":"Any"}
        ,{"id":"glycan", "display":"Glycan"}
        ,{"id":"protein", "display":"Protein"}
        ,{"id":"gene", "display":"Gene"}
        ,{"id":"organism", "display":"Organism"}
        ,{"id":"disease", "display":"Disease"}
        ,{"id":"pathway", "display":"Pathway"}
    ] 
    


    res_obj["aa_list"] = []
    for aathree in seen["amino_acid"]:
        if aathree not in aa_dict:
            continue
        res_obj["aa_list"].append(aa_dict[aathree])

    res_obj["glycosylation_evidence_type"] = [
        { "id": "all_sites", "display": "All sites"},
        { "id": "all_reported_sites_with_without_glycans","display": "All reported sites (with or without Glycans)"},
        { "id": "sites_reported_with_glycans", "display": "Sites reported with Glycans"},
        { "id": "sites_reported_without_glycans", "display": "Sites reported without Glycans"},
        { "id": "predicted_sites", "display": "Predicted sites"},
        { "id": "sites_detected_by_literature_mining", "display": "Sites detected by literature mining"}
    ]




    return res_obj

def get_aa_dict():

    aa_dict = {
        "Arg":{ "name":"Arginine - Arg - R", "key":"R" },
        "Lys":{ "name":"Lysine - Lys - K", "key":"K" },
        "Asp":{ "name":"Aspartic acid - Asp - D", "key":"D" },
        "Glu":{ "name":"Glutamic acid - Glu - E", "key":"E" },
        "Gln":{ "name":"Glutamine - Gln - Q", "key":"Q" },
        "Asn":{ "name":"Asparagine - Asn - N", "key":"N" },
        "His":{ "name":"Histidine - His - H", "key":"H" },
        "Ser":{ "name":"Serine - Ser - S", "key":"S" },
        "Thr":{ "name":"Threonine - Thr - T", "key":"T" },
        "Tyr":{ "name":"Tyrosine - Tyr - Y", "key":"Y" },
        "Cys":{ "name":"Cysteine - Cys - C", "key":"C" },
        "Trp":{ "name":"Tryptophan - Trp - W", "key":"W" },
        "Ala":{ "name":"Alanine - Ala - A", "key":"A" },
        "Ile":{ "name":"Isoleucine - Ile - I", "key":"I" },
        "Leu":{ "name":"Leucine - Leu - L", "key":"L" },
        "Met":{ "name":"Methionine - Met - M", "key":"M" },
        "Phe":{ "name":"Phenylalanine - Phe - F", "key":"F" },
        "Val":{ "name":"Valine - Val - V", "key":"V" },
        "Pro":{ "name":"Proline - Pro - P", "key":"P" },
        "Gly":{ "name":"Glycine - Gly - G", "key":"G" },
        "Sec":{ "name":"Selenocysteine - Sec - U", "key":"U"},
        "R":{ "name":"Arginine - Arg - R", "key":"R" },
        "K":{ "name":"Lysine - Lys - K", "key":"K" },
        "D":{ "name":"Aspartic acid - Asp - D", "key":"D" },
        "E":{ "name":"Glutamic acid - Glu - E", "key":"E" },
        "Q":{ "name":"Glutamine - Gln - Q", "key":"Q" },
        "N":{ "name":"Asparagine - Asn - N", "key":"N" },
        "H":{ "name":"Histidine - His - H", "key":"H" },
        "S":{ "name":"Serine - Ser - S", "key":"S" },
        "T":{ "name":"Threonine - Thr - T", "key":"T" },
        "Y":{ "name":"Tyrosine - Tyr - Y", "key":"Y" },
        "C":{ "name":"Cysteine - Cys - C", "key":"C" },
        "W":{ "name":"Tryptophan - Trp - W", "key":"W" },
        "A":{ "name":"Alanine - Ala - A", "key":"A" },
        "I":{ "name":"Isoleucine - Ile - I", "key":"I" },
        "L":{ "name":"Leucine - Leu - L", "key":"L" },
        "M":{ "name":"Methionine - Met - M", "key":"M" },
        "F":{ "name":"Phenylalanine - Phe - F", "key":"F" },
        "V":{ "name":"Valine - Val - V", "key":"V" },
        "P":{ "name":"Proline - Pro - P", "key":"P" },
        "G":{ "name":"Glycine - Gly - G", "key":"G" },
        "U":{ "name":"Selenocysteine - Sec - U", "key":"U" }
    }


    return aa_dict

###############################
def main():

        global config_obj
        global species_obj
        global species_list
        global aa_dict
        global PATTERN 


        DEBUG = False
        #DEBUG = True

        PATTERN = "*"
        if DEBUG:
            PATTERN = "P142*"

        config_obj = json.loads(open("../conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        
        aa_dict = get_aa_dict()

        cat_file = "generated/misc/simple_search_categories.json"
        cat_dict  = json.load(open(cat_file))
        

        species_obj = {}
        in_file = "generated/misc/species_info.csv"
        libgly.load_species_info(species_obj, in_file)

        species_list = []
        for k in species_obj:
            obj = species_obj[k]
            if obj["short_name"] not in species_list and obj["is_reference"] == "yes":
                species_list.append(obj["short_name"])

        log_file = "logs/make-searchinitdb.log"
        msg = "make-searchinitdb: stareted logging"
        csvutil.write_log_msg(log_file, msg, "w")


        out_obj = {}


        msg = "make-searchinitdb: ... started get_disease_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["disease"] = get_disease_search_init("c_disease")
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")
        #print (json.dumps(out_obj, indent=4))



        msg = "make-searchinitdb: ... started get_usecases_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["usecases"] = get_usecases_search_init()
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")

 
        msg = "make-searchinitdb: ... started get_biomarker_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["biomarker"] = get_biomarker_search_init("c_biomarker")
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")

            
        msg = "make-searchinitdb: ... started get_glycan_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["glycan"] = get_glycan_search_init("c_glycan")
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")

        msg = "make-searchinitdb: ... started get_idmapping_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["idmapping"] = get_idmapping_search_init()
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")

        msg = "make-searchinitdb: ... started get_site_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["site"] = get_site_search_init("c_site")
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")

        msg = "make-searchinitdb: ... started get_super_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["supersearch"] = get_super_search_init()
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")
            
        msg = "make-searchinitdb: ... started get_field_limits"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["fieldlimits"] = get_field_limits(config_obj["fieldlimitinfo"])
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")
            
        msg = "make-searchinitdb: ... started get_protein_search_init"
        csvutil.write_log_msg(log_file, msg, "a")
        out_obj["protein"] = get_protein_search_init("c_protein")
        msg = "make-searchinitdb: done!"
        csvutil.write_log_msg(log_file, msg, "a")


        for r_type in ["protein", "glycan", "biomarker"]:
            if r_type in out_obj and r_type in cat_dict:
                out_obj[r_type]["simple_search_category"] = cat_dict[r_type]
           
        out_file = "jsondb/searchinitdb/searchinitdb.json"
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(out_obj, indent=4)))
        msg = "make-searchinitdb: final created: %s searchinitdb objects" % (1)
        csvutil.write_log_msg(log_file, msg, "a")



if __name__ == '__main__':
	main()

