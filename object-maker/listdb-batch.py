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
import csvutil
import libgly
import subprocess


def extract_name(obj_list, name_type, resource):
   
    if obj_list == []:
        return ""
    
    name_list_dict = {"recommended":[], "synonym":[]}
    for obj in obj_list:
        if obj["resource"] == resource:
            name_list_dict[obj["type"]].append(obj["name"])


    
    if name_type == "recommended" and name_list_dict["recommended"] == []:
        name_list_dict["recommended"] += name_list_dict["synonym"]


    if name_type == "all":
        return "; ".join(name_list_dict["recommended"] + name_list_dict["synonym"])
    else:
        return "; ".join(name_list_dict[name_type])




def get_protein_record(obj):


    canon = obj["uniprot_canonical_ac"]
    uniprot_id = obj["uniprot_id"]
    protein_length = obj["sequence"]["length"]

    mass = -1
    if "mass" in obj:
        if "chemical_mass" in obj["mass"]:
            mass = obj["mass"]["chemical_mass"]
    

    protein_name, protein_names_uniprotkb, protein_names_refseq = "", "", ""
    if "protein_names" in obj:
        protein_name = extract_name(obj["protein_names"], "recommended", "UniProtKB")
        protein_names_uniprotkb = extract_name(obj["protein_names"], "all", "UniProtKB")
        protein_names_refseq = extract_name(obj["protein_names"], "all", "RefSeq")


    refseq_ac, refseq_name = "", ""
    if "refseq" in obj:
        refseq_ac = obj["refseq"]["ac"] if "ac" in obj["refseq"] else ""
        refseq_name = obj["refseq"]["name"] if "name" in obj["refseq"] else ""

    gene_name, gene_names_uniprotkb, gene_names_refseq = "", "", ""
    if "gene_names" in obj:
        gene_name = extract_name(obj["gene_names"], "recommended", "UniProtKB")
        gene_names_uniprotkb = extract_name(obj["gene_names"], "all", "UniProtKB")
        gene_names_refseq = extract_name(obj["gene_names"], "all", "RefSeq")



    species_map = json.loads(open("generated/misc/species_map.json", "r").read())
    organism_name, tax_id = "", 0
    if "species" in obj:
        if len(obj["species"]) > 0:
            organism_name = obj["species"][0]["glygen_name"] if "glygen_name" in obj["species"][0] else ""
            tax_id = obj["species"][0]["taxid"] if "taxid" in obj["species"][0] else 0
            tax_id_str = str(tax_id)
            if tax_id_str in species_map:
                organism_name = species_map[tax_id_str]["ref_tax_glygen_name"]

    disease_list = []
    if "disease" in obj:
        for o in obj["disease"]:
            disease_id = o["recommended_name"]["id"]
            disease_name = o["recommended_name"]["name"]
            resource = o["recommended_name"]["resource"]
            disease_list.append("%s (%s:%s)" % (disease_name, resource, disease_id))
    disease = "; ".join(disease_list)


    pathway_list = []
    if "pathway" in obj:
        for o in obj["pathway"]:
            pathway_id = o["id"]
            pathway_name = o["name"]
            pathway_list.append("%s (%s)" % (pathway_name, pathway_id))
    pathway = "; ".join(pathway_list)


    tmp_seen = {}
    if "biomarkers" in obj:
        for o in obj["biomarkers"]:
            for role in o["best_biomarker_role"]:
                tmp_seen[role] = True
    biomarker_type = "; ".join(list(tmp_seen.keys()))


    function_list = []
    if "function" in obj:
        for o in obj["function"]:
            ann = o["annotation"]
            badge_list= []
            for e in o["evidence"]:
                badge_list.append(e["database"])
            function_list.append("%s, [%s]" % (ann,", ".join(badge_list)))
    function = "; ".join(function_list)
    

    seen_publication_id = {}
    publication_count = 0
    if "publication" in obj:
        publication_count = len(obj["publication"])
        for o in obj["publication"]:
            for oo in o["reference"]:
                seen_publication_id[oo["id"]] = True

 
    seen_glycan = {} 
    predicted_glycosites, predicted_n_glycosites, predicted_o_glycosites = 0, 0, 0
    mined_glycosites, mined_n_glycosites, mined_o_glycosites = 0, 0, 0

    reported_n_glycosites, reported_o_glycosites = 0,0
    reported_n_glycosites_with_glycan, reported_o_glycosites_with_glycan = 0,0
    seen_fully_resolved = {}
    seen_site_id = {}
    missing_score = 10000
    if "glycosylation" in obj:
        for o in obj["glycosylation"]:
            start_pos = o["start_pos"] if "start_pos" in o else ""
            end_pos = o["end_pos"] if "end_pos" in o else ""
            glytoucan_ac = o["glytoucan_ac"]
            if glytoucan_ac != "":
                seen_glycan[glytoucan_ac] = True
                if glytoucan_ac in missing_score_dict:
                    if missing_score_dict[glytoucan_ac] < missing_score:
                        missing_score = missing_score_dict[glytoucan_ac]
                if glytoucan_ac in fully_resolved_dict:
                    seen_fully_resolved[glytoucan_ac] = True
            site_id = "%s.%s.%s" % (canon,start_pos, end_pos)
            site_type = o["type"].lower()
            if start_pos != "" and end_pos != "":
                seen_site_id[site_id] = True


            if o["site_category"] == "predicted":
                predicted_glycosites += 1
                if site_type == "n-linked":
                    predicted_n_glycosites += 1
                if site_type == "o-linked":
                    predicted_o_glycosites += 1
            elif site_type == "n-linked" and o["site_category"] == "reported":
                reported_n_glycosites += 1
            elif site_type == "o-linked" and o["site_category"] == "reported":
                reported_o_glycosites += 1
            elif site_type == "n-linked" and o["site_category"] == "reported_with_glycan":
                reported_n_glycosites_with_glycan += 1
            elif site_type == "o-linked" and o["site_category"] == "reported_with_glycan":
                reported_o_glycosites_with_glycan += 1
            
            if "automatic_literature_mining" in o["site_category_dict"]:
                mined_glycosites += 1
                if site_type == "n-linked":
                    mined_n_glycosites += 1
                if site_type == "o-linked":
                    mined_o_glycosites += 1   
 
    reported_phosphosites = 0
    if "phosphorylation" in obj:
        reported_phosphosites = len(obj["phosphorylation"])
    
    repoted_snv = 0
    if "snv" in obj:
        reported_snv = len(obj["snv"])
    
    reported_mutagensis = 0
    if "mutagenesis" in obj:
        reported_mutagensis = len(obj["mutagenesis"])
    
    reported_glycation = 0
    if "glycation" in obj:
        reported_glycation = len(obj["glycation"])
   
    seen_binding_glycan = {} 
    reported_interactions = 0
    if "interactions" in obj:
        reported_interactions = len(obj["interactions"])
        for o in obj["interactions"]:
            seen_binding_glycan[o["interactor_id"]] = True

    seen_ec_number = {}
    if "reaction_enzymes" in obj:
        for o in obj["reaction_enzymes"]:
            if "activity" in o:
                for oo in o["activity"]:
                    seen_ec_number[oo["ec_number"]] = True


    seen_synthesized_glycans = {}
    if "synthesized_glycans" in obj:
        for o in obj["synthesized_glycans"]:
            seen_synthesized_glycans[o["glytoucan_ac"]] = True
    seen_go_ann = {}
    if "go_annotation" in obj:
        if "categories" in obj["go_annotation"]:
            for o in obj["go_annotation"]["categories"]:
                cat = o["name"]
                if "go_terms" in o:
                    for oo in o["go_terms"]:
                        go_term = ""
                        if "id" in oo and "name" in oo:
                            go_term = oo["id"] + " " + oo["name"]
                            if cat not in seen_go_ann:
                                seen_go_ann[cat] = []
                            if len(seen_go_ann[cat]) < 5:
                                seen_go_ann[cat].append(go_term)


    total_reported_n_glycosites = reported_n_glycosites + reported_n_glycosites_with_glycan
    total_reported_o_glycosites = reported_o_glycosites + reported_o_glycosites_with_glycan
    total_n_glycosites = total_reported_n_glycosites + predicted_n_glycosites
    total_o_glycosites = total_reported_o_glycosites + predicted_o_glycosites

    go_f, go_p, go_c = "", "", ""
    cat = "Molecular Function"
    if cat in seen_go_ann:
        go_f = ";".join(seen_go_ann[cat])
    cat = "Cellular Component"
    if cat in seen_go_ann:
        go_c = ";".join(seen_go_ann[cat])
    cat = "Biological Process"
    if cat in seen_go_ann:
        go_p = ";".join(seen_go_ann[cat])

    





    glycan_count = len(seen_glycan.keys())
    disease_count = len(disease_list)
    publication_id_list = ";".join(list(seen_publication_id.keys()))
    glyco_site_count = len(seen_site_id.keys())
    binding_glycan_count = len(seen_binding_glycan.keys())
    ec_number_list = ";".join(list(seen_ec_number.keys()))
    #synthesized_glycan_list = ";".join(list(seen_synthesized_glycans.keys()))
    synthesized_glycan_count = len(seen_synthesized_glycans.keys())

    uniprot_ac = canon.split("-")[0]
    secondary_ac_list = secondary_ac_dict[uniprot_ac] if uniprot_ac in secondary_ac_dict else ""


    o = {
        "uniprot_canonical_ac":canon
        ,"uniprot_id":uniprot_id
        ,"mass": mass
        ,"length":protein_length
        ,"protein_name":protein_name
        ,"gene_name":gene_name
        ,"organism":organism_name
        ,"refseq_name": refseq_name
        ,"refseq_ac": refseq_ac
        ,"gene_names_uniprotkb":gene_names_uniprotkb
        ,"gene_names_refseq":gene_names_refseq
        ,"protein_names_uniprotkb":protein_names_uniprotkb
        ,"protein_names_refseq":protein_names_refseq
        ,"tax_id":tax_id
        ,"disease":disease
        ,"pathway":pathway
        ,"pathway_count":len(pathway_list)
        ,"publication_count":publication_count
        ,"publication_id_list":publication_id_list
        ,"function":function
        ,"predicted_glycosites":predicted_glycosites
        ,"predicted_n_glycosites":predicted_n_glycosites
        ,"predicted_o_glycosites":predicted_o_glycosites
        ,"mined_glycosites":mined_glycosites
        ,"mined_n_glycosites":mined_n_glycosites
        ,"mined_o_glycosites":mined_o_glycosites
        ,"reported_n_glycosites":reported_n_glycosites
        ,"reported_o_glycosites":reported_o_glycosites
        ,"reported_n_glycosites_with_glycan":reported_n_glycosites_with_glycan
        ,"reported_o_glycosites_with_glycan":reported_o_glycosites_with_glycan
        ,"total_reported_n_glycosites":total_reported_n_glycosites
        ,"total_reported_o_glycosites":total_reported_o_glycosites
        ,"total_n_glycosites": total_n_glycosites
        ,"total_o_glycosites": total_o_glycosites
        ,"reported_fully_resolved_glycans":len(list(seen_fully_resolved.keys()))
        ,"reported_phosphosites":reported_phosphosites
        ,"reported_mutagensis":reported_mutagensis
        ,"reported_glycation":reported_glycation
        ,"reported_snv":reported_snv
        ,"reported_interactions":reported_interactions
        ,"missing_score":missing_score
        ,"biomarker_type":biomarker_type
        ,"glycan_count":glycan_count
        ,"disease_count":disease_count
        ,"glyco_site_count":glyco_site_count 
        ,"binding_glycan_count":binding_glycan_count
        ,"secondary_ac_list":secondary_ac_list
        ,"ec_number_list":ec_number_list
        ,"go_molecular_function":go_f
        ,"go_cellular_component":go_c 
        ,"go_biological_process":go_p   
        ,"synthesized_glycan_count":synthesized_glycan_count
    }



    return o


def get_motif_record(doc):
   

    seen = {"protein":{}, "enzyme":{}}
    for o in doc["glycans"]:
        glytoucan_ac = o["glytoucan_ac"]
        json_file = "jsondb/glycandb/%s.json" % (glytoucan_ac)
        glycan_doc = json.loads(open(json_file,"r").read())
        for oo in glycan_doc["glycoprotein"]:
            uniprot_canonical_ac = oo["uniprot_canonical_ac"]
            seen["protein"][uniprot_canonical_ac] = True
        #for oo in glycan_doc["enzyme"]:
        #    uniprot_canonical_ac = oo["uniprot_canonical_ac"]
        #    seen["protein"][uniprot_canonical_ac] = True
        #    seen["enzyme"][uniprot_canonical_ac] = True


    
    tmp_seen = {}
    if "biomarkers" in doc:
        for o in doc["biomarkers"]:
            for role in o["best_biomarker_role"]:
                tmp_seen[role] = True
    biomarker_type = "; ".join(list(tmp_seen.keys()))


    protein_count =  len(list(seen["protein"].keys()))
    enzyme_count = len(list(seen["enzyme"].keys()))
    obj =  {"motif_ac": doc["motif_ac"], "motif_name":"",
        "glytoucan_ac": doc["glytoucan_ac"],
        "glycan_count": len(doc["glycans"]), 
        "protein_count": protein_count,
        "enzyme_count":enzyme_count,
        "biomarker_type":biomarker_type,
        "synonyms":[], 
        "pmidlist":[]
    }
    if len(doc["names"]) > 0:
        obj["motif_name"] = doc["names"][0]["name"]
    for o in doc["publication"]:
        for oo in o["reference"]:
            if oo["type"] == "PUBMED":
                obj["pmidlist"].append(oo["id"])
    
    obj["synonyms"] = doc["synonym"]


    return obj


def get_glycan_record(doc):

    glytoucan_ac = doc["glytoucan_ac"] 
    seen_connector = {"uniprot_canonical_ac":{}, "enzyme":{}, "site":{}, "species":{}, 
            "motif":{}, "interactions":{}, "sequence_details":{}}


    if "glycoprotein" in doc:
        xobj_list = doc["glycoprotein"]
        file_list = glob.glob("jsondb/batchdb/glycan.%s.*.json" % (glytoucan_ac))
        for in_file in file_list:
            batch_doc = json.loads(open(in_file, "r").read())
            if "glycoprotein" in batch_doc["sections"]:
                for xobj in batch_doc["sections"]["glycoprotein"]:
                    xobj_list.append(xobj) 
        for xobj in xobj_list:
            canon = xobj["uniprot_canonical_ac"]
            combo_id = canon
            seen_connector["uniprot_canonical_ac"][combo_id] = True
            if "position" in xobj:
                start_pos, end_pos = xobj["position"], xobj["position"]
                combo_id = "%s.%s.%s" % (canon,start_pos, end_pos)
                seen_connector["site"][combo_id] = True
   
    
    tmp_seen = {}
    if "biomarkers" in doc:
        for o in doc["biomarkers"]:
            for role in o["best_biomarker_role"]:
                tmp_seen[role] = True

    biomarker_type = "; ".join(list(tmp_seen.keys()))


    if "glycan_type" in doc:
        if doc["glycan_type"] not in ["", "Saccharide"]:
            seen_connector["sequence_details"][doc["glycan_type"]] = True

    seen_name = {}
    if "names" in doc:
        for o in doc["names"]:
            name = "%s[%s]" % (o["name"], o["domain"])
            seen_name[name] = True
    if "fully_determined" in doc:
        if doc["fully_determined"].lower() == "yes":
            seen_connector["sequence_details"]["Fully defined"] = True
    
    list_one = ["Fully defined", "Topology", "Composition", "BaseComposition"]
    list_two = list(seen_connector["sequence_details"].keys())
    list_three = list(set(list_one).intersection(set(list_two)))
    
    if list_three == []:
        seen_connector["sequence_details"]["Incomplete"] = True



    if "species" in doc:
        species_map = json.loads(open("generated/misc/species_map.json", "r").read())
        for xobj in doc["species"]:
            name = xobj["glygen_name"]
            tax_id = xobj["taxid"] if "taxid" in xobj else 0
            tax_id_str = str(tax_id)
            if tax_id_str in species_map:
                name = species_map[tax_id_str]["ref_tax_glygen_name"]
            seen_connector["species"][name] = True



    if "motifs" in doc:
        for xobj in doc["motifs"]:
            seen_connector["motif"][xobj["id"]] = True

    seen_interactor_id = {}
    if "interactions" in doc:
        for xobj in doc["interactions"]:
            seen_connector["interactions"][xobj["interactor_id"]] = True
            
    seen_pubchem_id, seen_chebi_id = {}, {}
    if "crossref" in doc:
        for o in doc["crossref"]:
            if o["database"].find("PubChem") != -1:
                seen_pubchem_id[o["id"]] = True
            if o["database"].find("ChEBI") != -1:
                seen_chebi_id[o["id"]] = True

    if "enzyme" in doc:
        for xobj in doc["enzyme"]:
            if "uniprot_canonical_ac" in xobj:
                seen_connector["enzyme"][xobj["uniprot_canonical_ac"]] = True


    seen_comp_id, seen_base_comp_id = {}, {}
    if "subsumption" in doc:
        for o in doc["subsumption"]:
            if o["relationship"] == "composition":
                seen_comp_id[o["related_accession"]] = True
            if o["relationship"] == "basecomposition":
                seen_base_comp_id[o["related_accession"]] = True
                




    missing_score = doc["missing_score"] 
    keywords = ";".join(doc["keywords"]) if "keywords" in doc else ""
    mass = doc["mass"] if "mass" in doc else -1
    mass_pme = doc["mass_pme"] if "mass_pme" in doc else -1
    number_monosaccharides = -1
    if "number_monosaccharides" in doc:
        number_monosaccharides = doc["number_monosaccharides"]
    seq_dict = {}
    
    image_url = "https://api.glygen.org/glycan/image/%s" % (glytoucan_ac)
    record_obj = {
        "glytoucan_ac":glytoucan_ac
        ,"image_url":image_url
        ,"names":";".join(list(seen_name.keys()))
        ,"mass": mass 
        ,"mass_pme": mass_pme
        ,"number_monosaccharides": number_monosaccharides
        ,"number_proteins": len(list(seen_connector["uniprot_canonical_ac"].keys()))
        ,"number_enzymes": len(list(seen_connector["enzyme"].keys()))
        ,"number_species": len(list(seen_connector["species"].keys()))
        ,"interactions":";".join(list(seen_connector["interactions"].keys()))
        ,"interactions_count":len(seen_connector["interactions"].keys())
        ,"species_list":";".join(list(seen_connector["species"].keys()))
        ,"sequence_details":";".join(list(seen_connector["sequence_details"].keys()))
        ,"pubchem_id":";".join(list(seen_pubchem_id.keys())) 
        ,"chebi_id":";".join(list(seen_chebi_id.keys()))
        ,"composition_id":";".join(list(seen_comp_id.keys()))
        ,"base_composition_id":";".join(list(seen_base_comp_id.keys()))
        ,"keywords":keywords
        ,"hit_score":0.0
        ,"missing_score":missing_score
        ,"biomarker_type":biomarker_type
    }
    
    for k in ["iupac", "glycoct", "glycam", "wurcs", "inchi", "inchi_key", "byonic", "gwb"]:
        record_obj[k] = ""
        if k in doc:
            record_obj[k] =  doc[k]
            if k == "inchi_key":
                record_obj[k] = ""
                if "key" in doc[k]:
                    record_obj[k] = doc[k]["key"]
    p_dict = {
        "enzyme":"uniprot_canonical_ac", 
        "glycoprotein":"uniprot_canonical_ac",
        "motifs":"id",
        "composition":"residue",
        "publication":"",
        "classification":""
    }
    for k in p_dict:
        tmp_dict = {}
        if k in doc:
            for oo in doc[k]:
                if k == "composition":
                    if oo["count"] > 0:
                            val = "%s (%s)" % (oo["residue"], oo["count"])
                            tmp_dict[val] = True
                elif k == "publication":
                    if "reference" in oo:
                        val = oo["reference"][0]["id"]
                        tmp_dict[val] = True
                elif k == "classification":
                    val = "%s" %(oo["type"]["name"])
                    tmp_dict[val] = True
                    val = "%s/%s" %(oo["type"]["name"], oo["subtype"]["name"])
                    tmp_dict[val] = True
                else:
                    val = oo[p_dict[k]]
                    tmp_dict[val] = True
        record_obj[k] = "; ".join(list(tmp_dict.keys()))
        if k in ["glycoprotein", "publication"]:
            kk = k + "_count"
            record_obj[kk] = len(tmp_dict.keys())

    for k in ["mass", "mass_pme", "number_monosaccharides"]:
        if record_obj[k] == -1:
            record_obj.pop(k)
   

 
    return record_obj



def get_biomarker_record(doc):

    record_obj = {
        "biomarker_canonical_id":doc["biomarker_canonical_id"],
        "biomarker_id":doc["biomarker_id"]
    }
    f_list = ["biomarker", "assessed_biomarker_entity_id", "assessed_entity_type"]
    f_list += ["assessed_biomarker_entity", "specimen", "best_biomarker_role"]
    f_list += ["condition"]
    seen = {}
    for f in f_list:
        seen[f] = {}


    seen_cmp_evdn = {}
    cmp_count = 0
    for obj in doc["biomarker_component"]:
        cmp_count += 1
        for f in ["biomarker", "assessed_biomarker_entity_id", "assessed_entity_type"]:
            val = obj[f]
            seen[f][val] = True
        f = "assessed_biomarker_entity"
        val = obj[f]["recommended_name"]
        seen[f][val] = True 
        f = "specimen"
        for o in obj[f]:
            val = "%s (%s)" % (o["name"], o["id"])
            seen[f][val] = True
        f = "evidence_source"
        for o in obj[f]:
            if o["database"] in ["PubMed", "DOI"]:
                seen_cmp_evdn[o["id"]] = True

    f = "best_biomarker_role"
    for obj in doc[f]:
        val = obj["role"]
        seen[f][val] = True

    f = "condition"
    if f in doc:
        val = "%s (%s)"%(doc[f]["recommended_name"]["name"], doc[f]["recommended_name"]["id"])
        seen[f][val] = True
        #for o in doc[f]["synonyms"]:
        #    val = "%s (%s)" % (o["name"], o["id"])
        #    seen[f][val] = True

    seen_glbl_evdn = {}
    f = "evidence_source"
    for o in doc[f]:
        if o["database"] in ["PubMed", "DOI"]:
            seen_glbl_evdn[o["id"]] = True


    for f in f_list:
        if f not in seen:
            seen[f][""] = True
        record_obj[f] = "; ".join(seen[f].keys())

    record_obj["component_count"] = cmp_count
    record_obj["component_evidence_count"] = len(list(seen_cmp_evdn.keys()))
    record_obj["global_evidence_count"] = len(list(seen_glbl_evdn.keys()))


    #print (json.dumps(record_obj, indent=4))
    return record_obj


def get_disease_record(doc):

    url_map = {
        "doid":"http://disease-ontology.org/term/DOID:%s",
        "mondo":"https://monarchinitiative.org/disease/MONDO:%s",
        "mim":"http://www.omim.org/entry/%s"
    }

    record_obj = {
        "record_id":doc["record_id"],
        "disease_id":doc["disease_id"],
        "recommended_name":doc["recommended_name"]["name"],    
        "record_url":""
    }

    prefix, id_value = doc["record_id"].split(".")[0], doc["record_id"].split(".")[1]
    if prefix in url_map:
        record_obj["record_url"] = url_map[prefix] % (id_value)   

    seen = {}
    seen["snv_count"] = 0
    seen["protein_count"] = len(doc["proteins"]) 
    seen["glycan_count"] = len(doc["glycans"])
    seen["biomarker_count"] = len(doc["biomarkers"])
     
    

    seen_species = {}
    for obj in doc["proteins"] + doc["glycans"]:
        for o in obj["species"]:
            name = o["glygen_name"]
            seen_species[name] = True
    seen["species_list"] = ";".join(list(seen_species.keys()))

    seen_biomarker_role = {}
    if len(doc["biomarkers"]) > 0:
        for obj in doc["biomarkers"]:
            for o in obj["best_biomarker_role"]:
                role = o["role"]
                seen_biomarker_role[role] = True

    seen["biomarker_type"] = ";".join(list(seen_biomarker_role.keys()))


    for f in seen:
        record_obj[f] = seen[f]

    return record_obj
    






def get_site_record(doc):

    json_file = "jsondb/proteindb/%s.json" % (doc["uniprot_canonical_ac"])
    if os.path.isfile(json_file) == False:
        return {}

    protein_doc = json.loads(open(json_file,"r").read())
    protein_name = extract_name(protein_doc["protein_names"], "recommended", "UniProtKB")
    uniprot_id = protein_doc["uniprot_id"]
    refseq_ac, refseq_name = "", ""
    if "refseq" in protein_doc:
        refseq_ac = protein_doc["refseq"]["ac"] if "ac" in protein_doc["refseq"] else ""
        refseq_name = protein_doc["refseq"]["name"] if "name" in protein_doc["refseq"] else ""
   
    gene_name, gene_names_uniprotkb, gene_names_refseq = "", "", ""
    if "gene_names" in protein_doc:
        gene_name = extract_name(protein_doc["gene_names"], "recommended", "UniProtKB")
        gene_names_uniprotkb = extract_name(protein_doc["gene_names"], "all", "UniProtKB")
        gene_names_refseq = extract_name(protein_doc["gene_names"], "all", "RefSeq")

    protein_length = protein_doc["sequence"]["length"]
 

 
    site_seq = doc["site_seq"] if "site_seq" in doc else ""
    up_seq = doc["up_seq"] if "up_seq" in doc else ""
    down_seq = doc["down_seq"] if "down_seq" in doc else ""



    record_obj = {
        "uniprot_canonical_ac":doc["uniprot_canonical_ac"]
        ,"protein_name":protein_name
        ,"start_pos":doc["start_pos"]
        ,"end_pos":doc["end_pos"]
        ,"site_seq":site_seq
        ,"up_seq":up_seq
        ,"down_seq":down_seq
        ,"uniprot_id":uniprot_id
        ,"refseq_ac":refseq_ac
        ,"refseq_name":refseq_name
        ,"gene_name":gene_name
        ,"length":protein_length
    }
    sec_list = config_obj["sitesections"]
    missing_score = 10000 
    seen = {"residue":{}, "glycosylation_type":{}, "snv_type":{}}
    n_glycosylation_count, o_glycosylation_count = 0, 0
    reported_flag, mined_flag, predicted_flag = False, False, False
    seen_glycan = {}
    for sec in sec_list:
        if sec in doc:
            record_obj[sec] = "yes" if len(doc[sec]) > 0 else "no"
            record_obj[sec+"_count"] = len(doc[sec]) 
            for obj in doc[sec]:
                if "residue" in obj:
                    seen["residue"][obj["residue"]] = True
                if sec == "glycosylation":
                    gtype = obj["type"].lower() if "type" in obj else ""
                    n_glycosylation_count += 1 if gtype == "n-linked" else 0
                    o_glycosylation_count += 1 if gtype == "o-linked" else 0
                    gsubtype = obj["subtype"].lower() if "subtype" in obj else ""
                    if gtype != "" and gsubtype == "":
                        gsubtype = "other"
                    value = "%s|%s" % (gtype, gsubtype) if gtype != "" else ""
                    seen["glycosylation_type"][value] = True
                    glytoucan_ac = obj["glytoucan_ac"]
                    if glytoucan_ac != "":
                        seen_glycan[glytoucan_ac] = True
                    if glytoucan_ac in missing_score_dict:
                        if missing_score_dict[glytoucan_ac] < missing_score:
                            missing_score = missing_score_dict[glytoucan_ac]
                    
                    cat_text = ";".join(list(obj["site_category_dict"].keys()))
                    if cat_text.find("reported") != -1:
                        reported_flag = True
                    if cat_text.find("automatic_literature_mining") != -1: 
                        mined_flag = True
                    if cat_text.find("predicted") != -1: 
                        predicted_flag = True
 
                if sec == "snv":
                    for value in obj["keywords"]:
                        if value.lower() in ["germline", "somatic"]:
                            seen["snv_type"][value.lower()] = True

        else:
            record_obj[sec] = ""
  
    record_obj["residue"] = list(seen["residue"].keys())[0] if seen["residue"] != {} else ""
    record_obj["glycosylation_type"] = ";".join(list(seen["glycosylation_type"].keys()))
    record_obj["snv_type"] = ";".join(list(seen["snv_type"].keys()))
    record_obj["missing_score"] = missing_score

    record_obj["n_glycosylation_count"] = n_glycosylation_count
    record_obj["o_glycosylation_count"] = o_glycosylation_count
    record_obj["glycan_count"] = len(seen_glycan.keys())
    record_obj["reported_glycosite_flag"] = reported_flag
    record_obj["mined_glycosite_flag"] = mined_flag
    record_obj["predicted_glycosite_flag"] = predicted_flag
 
    record_obj["tax_id"], record_obj["organism"] = "", 0
    if "species" in doc:
        if len(doc["species"]) > 0:
            record_obj["organism"] = doc["species"][0]["glygen_name"] if "glygen_name" in doc["species"][0] else ""
            record_obj["tax_id"] = doc["species"][0]["taxid"] if "taxid" in doc["species"][0] else 0



    return record_obj









def main():

    
    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " )
    parser.add_option("-s","--start",action="store",dest="start",help="")
    parser.add_option("-e","--end",action="store",dest="end",help="")

    (options,args) = parser.parse_args()
    for file in ([options.start, options.end]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    start = int(options.start)
    end = int(options.end)






    global config_obj
    global path_obj
    global species_list
    global map_dict
    global data_dir
    global main_dict
    global missing_score_dict
    global secondary_ac_dict
    global fully_resolved_dict


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"

    #DEBUG = False
    DEBUG = True

    pattern_dict = {
        "protein":"*", "glycan":"*", "site":"*", "motif":"*", "biomarker":"*",
        "disease":"*"
    }
    if DEBUG:
        pattern_dict = {"disease":"*"}
        #pattern_dict = {"protein":"P04637*"}
        #pattern_dict = {"glycan":"*", "protein":"*"}
        #pattern_dict = {"glycan":"*"}
           


    file_obj_list = []
    for record_type in sorted(list(pattern_dict.keys())):
        #glob_str = "temp/%s.json" % (pattern_dict[record_type])
        glob_str = "jsondb/%sdb/%s.json" % (record_type, pattern_dict[record_type])
        file_list = glob.glob(glob_str)
        for json_file in sorted(file_list):
            file_obj_list.append({"file":json_file, "record_type":record_type})
   
    #print (len(file_obj_list))
    #exit()


    end = end if end < len(file_obj_list) else len(file_obj_list)
    
    log_file = "logs/make-listdb.%s-%s.log" % (start, end)
    msg = "make-listdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")


    species_list = []
    species_obj = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)
    for tax_id in species_obj:
        if species_obj[tax_id]["is_reference"] == "yes":
            name = species_obj[tax_id]["glygen_name"]
            if name not in species_list:
                species_list.append(name)
    species_list = sorted(species_list)

    # load some dicts
    missing_score_dict = {}
    secondary_ac_dict = {}
    fully_resolved_dict = {}
    if DEBUG == False:
        msg = "make-listdb: loading secondary_ac_dict"
        csvutil.write_log_msg(log_file, msg, "a")
        file_list = glob.glob("downloads/ebi/current/accession-history-*.tsv")
        for in_file in file_list:
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, "\t")
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                ac = row[-1]
                secondary_ac_list = row[0].replace(" ", "").split(",")
                secondary_ac_list = ";".join(secondary_ac_list)
                secondary_ac_dict[ac] = secondary_ac_list

        msg = "make-listdb: loading fully_resolved_dict and missing_score_dict"
        csvutil.write_log_msg(log_file, msg, "a")
        file_list = glob.glob("jsondb/glycandb/*.json")
        file_list += glob.glob("jsondb/jumbodb/glycandb/*.json")
        for json_file in file_list:
            glytoucan_ac = json_file.split("/")[-1].split(".")[0]
            doc = json.loads(open(json_file,"r").read())
            missing_score_dict[glytoucan_ac] = doc["missing_score"]
            if doc["missing_score"] == 0:
                fully_resolved_dict[glytoucan_ac] = True


    record_count_dict = {}
    for obj in file_obj_list[start-1:end]:
        record_type = obj["record_type"]
        json_file = obj["file"]
        if record_type not in record_count_dict:
            record_count_dict[record_type] = 0
        r = record_count_dict[record_type]
        if r > 0 and r%1000 == 0:
            msg = "make-listdb: created %s %s listdb records" % (record_type, r)
            csvutil.write_log_msg(log_file, msg, "a")
        doc = json.loads(open(json_file,"r").read())
        record_obj = {}
        if record_type == "protein":
            record_obj = get_protein_record(doc)
            record_id = "%s" % (doc["uniprot_canonical_ac"])
        elif record_type == "glycan":
            record_obj = get_glycan_record(doc)
            record_id = "%s" % (doc["glytoucan_ac"])
        elif record_type == "motif":
            record_obj = get_motif_record(doc)
            record_id = "%s" % (doc["motif_ac"])
        elif record_type == "site":
            record_obj = get_site_record(doc)
            record_id = "%s" % (doc["id"])
        elif record_type == "biomarker":
            record_obj = get_biomarker_record(doc)
            record_id = "%s" % (doc["biomarker_id"])
        elif record_type == "disease":
            record_obj = get_disease_record(doc)
            record_id = "%s" % (doc["record_id"])

        if record_obj == {}:
            continue
        record_obj["record_id"] = record_id
        record_obj["record_type"] = record_type
        out_file = "jsondb/listdb/%s.json" % (record_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(record_obj, indent=4)))
        #if DEBUG:
        #    print ("created %s" % (out_file))
        if DEBUG and r%100 == 0:
            msg = "make-listdb: final created: %s %slist records" % (r, record_type)
            #print (msg)
        record_count_dict[record_type] += 1


    for record_type in record_count_dict:
        r = record_count_dict[record_type] 
        msg = "make-listdb: final created: %s %slist records" % (r, record_type)
        csvutil.write_log_msg(log_file, msg, "a")




if __name__ == '__main__':
    main()

