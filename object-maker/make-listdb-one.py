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


    organism_name, tax_id = "", 0
    if "species" in obj:
        if len(obj["species"]) > 0:
            organism_name = obj["species"][0]["name"] if "name" in obj["species"][0] else ""
            tax_id = obj["species"][0]["taxid"] if "taxid" in obj["species"][0] else 0

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


    function_list = []
    if "function" in obj:
        for o in obj["function"]:
            ann = o["annotation"]
            badge_list= []
            for e in o["evidence"]:
                badge_list.append(e["database"])
            function_list.append("%s, [%s]" % (ann,", ".join(badge_list)))
    function = "; ".join(function_list)
    
    publication_count = 0
    if "publication" in obj:
        publication_count = len(obj["publication"])
  
    predicted_glycosites =  0
    reported_n_glycosites, reported_o_glycosites = 0,0
    reported_n_glycosites_with_glycan, reported_o_glycosites_with_glycan = 0,0
    seen_fully_resolved = {}
    missing_score = 10000
    if "glycosylation" in obj:
        for o in obj["glycosylation"]:
            aa_pos = o["position"] if "position" in o else ""
            glytoucan_ac = o["glytoucan_ac"]
            if glytoucan_ac != "":
                if glytoucan_ac in missing_score_dict:
                    if missing_score_dict[glytoucan_ac] < missing_score:
                        missing_score = missing_score_dict[glytoucan_ac]
                json_file = "jsondb/glycandb/%s.json" % (glytoucan_ac)
                glycan_doc = json.loads(open(json_file,"r").read())
                if glycan_doc["missing_score"] == 0:
                    seen_fully_resolved[glytoucan_ac] = True
            
            site_id = "%s.%s.%s" % (canon,aa_pos, aa_pos)
            site_type = o["type"].lower()

            if o["site_category"] == "predicted":
                predicted_glycosites += 1
            elif site_type == "n-linked" and o["site_category"] == "reported":
                reported_n_glycosites += 1
            elif site_type == "o-linked" and o["site_category"] == "reported":
                reported_o_glycosites += 1
            elif site_type == "n-linked" and o["site_category"] == "reported_with_glycan":
                reported_n_glycosites_with_glycan += 1
            elif site_type == "o-linked" and o["site_category"] == "reported_with_glycan":
                reported_o_glycosites_with_glycan += 1
    
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
    
    reported_interactions = 0
    if "interactions" in obj:
        reported_interactions = len(obj["interactions"])

    o = {
        "uniprot_canonical_ac":canon
        ,"mass": mass
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
        ,"publication_count":publication_count
        ,"function":function
        ,"predicted_glycosites":predicted_glycosites
        ,"reported_n_glycosites":reported_n_glycosites
        ,"reported_o_glycosites":reported_o_glycosites
        ,"reported_n_glycosites_with_glycan":reported_n_glycosites_with_glycan
        ,"reported_o_glycosites_with_glycan":reported_o_glycosites_with_glycan
        ,"reported_fully_resolved_glycans":len(seen_fully_resolved.keys())
        ,"reported_phosphosites":reported_phosphosites
        ,"reported_mutagensis":reported_mutagensis
        ,"reported_glycation":reported_glycation
        ,"reported_snv":reported_snv
        ,"reported_interactions":reported_interactions
        ,"missing_score":missing_score
    }

    return o


def get_motif_record(doc):
   

    seen = {}
    for o in doc["glycans"]:
        glytoucan_ac = o["glytoucan_ac"]
        json_file = "jsondb/glycandb/%s.json" % (glytoucan_ac)
        glycan_doc = json.loads(open(json_file,"r").read())
        for oo in glycan_doc["glycoprotein"]:
            uniprot_canonical_ac = oo["uniprot_canonical_ac"]
            seen[uniprot_canonical_ac] = True
    
    protein_count =  len(seen.keys())
    obj =  {"motif_ac": doc["motif_ac"], "motif_name":"",
        "glytoucan_ac": doc["glytoucan_ac"],
        "glycan_count": len(doc["glycans"]), 
        "protein_count": protein_count,
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
        for xobj in doc["glycoprotein"]:
            canon = xobj["uniprot_canonical_ac"]
            combo_id = canon
            seen_connector["uniprot_canonical_ac"][combo_id] = True
            if "position" in xobj:
                start_pos, end_pos = xobj["position"], xobj["position"]
                combo_id = "%s.%s.%s" % (canon,start_pos, end_pos)
                seen_connector["site"][combo_id] = True
   

    if "glycan_type" in doc:
        if doc["glycan_type"] not in ["", "Saccharide"]:
            seen_connector["sequence_details"][doc["glycan_type"]] = True

    
    if "fully_determined" in doc:
        if doc["fully_determined"].lower() == "yes":
            seen_connector["sequence_details"]["Fully defined"] = True
    
    list_one = ["Fully defined", "Topology", "Composition", "BaseComposition"]
    list_two = seen_connector["sequence_details"].keys()
    list_three = list(set(list_one).intersection(set(list_two)))
    
    if list_three == []:
        seen_connector["sequence_details"]["Incomplete"] = True



    if "species" in doc:
        for xobj in doc["species"]:
            seen_connector["species"][xobj["name"]] = True

    if "motifs" in doc:
        for xobj in doc["motifs"]:
            seen_connector["motif"][xobj["id"]] = True

    if "interactions" in doc:
        for xobj in doc["interactions"]:
            seen_connector["interactions"][xobj["interactor_id"]] = True



    if "enzyme" in doc:
        for xobj in doc["enzyme"]:
            if "uniprot_canonical_ac" in xobj:
                seen_connector["enzyme"][xobj["uniprot_canonical_ac"]] = True

    missing_score = doc["missing_score"] 
    keywords = ";".join(doc["keywords"]) if "keywords" in doc else ""
    mass = doc["mass"] if "mass" in doc else -1
    mass_pme = doc["mass_pme"] if "mass_pme" in doc else -1
    number_monosaccharides = -1
    if "number_monosaccharides" in doc:
        number_monosaccharides = doc["number_monosaccharides"]
    seq_dict = {}
    record_obj = {
        "glytoucan_ac":glytoucan_ac
        ,"mass": mass 
        ,"mass_pme": mass_pme
        ,"number_monosaccharides": number_monosaccharides
        ,"number_proteins": len(seen_connector["uniprot_canonical_ac"].keys())
        ,"number_enzymes": len(seen_connector["enzyme"].keys())
        ,"number_species": len(seen_connector["species"].keys())
        ,"interactions":";".join(seen_connector["interactions"].keys())
        ,"species_list":";".join(seen_connector["species"].keys())
        ,"sequence_details":";".join(seen_connector["sequence_details"].keys())
        ,"keywords":keywords
        ,"hit_score":0.0
        ,"missing_score":missing_score
    }

    for k in ["iupac", "glycoct", "glycam", "wurcs", "inchi", "inchi_key", "byonic", "gwb"]:
        record_obj[k] = ""
        if k in doc:
            record_obj[k] =  doc[k]
            if k == "inchi_key":
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
        tmp_list = []
        if k in doc:
            for oo in doc[k]:
                if k == "composition":
                    if oo["count"] > 0:
                            val = "%s (%s)" % (oo["residue"], oo["count"])
                            tmp_list.append(val)
                elif k == "publication":
                    if "reference" in oo:
                        val = oo["reference"][0]["id"]
                        tmp_list.append(val)
                elif k == "classification":
                    val = "%s/%s" %(oo["type"]["name"], oo["subtype"]["name"])
                    tmp_list.append(val)
                else:
                    val = oo[p_dict[k]]
                    tmp_list.append(val)
        record_obj[k] = "; ".join(tmp_list)
        
    for k in ["mass", "mass_pme", "number_monosaccharides"]:
        if record_obj[k] == -1:
            record_obj.pop(k)
    
    return record_obj







def get_site_record(doc):

    json_file = "jsondb/proteindb/%s.json" % (doc["uniprot_canonical_ac"])
    protein_doc = json.loads(open(json_file,"r").read())
    protein_name = extract_name(protein_doc["protein_names"], "recommended", "UniProtKB")

    record_obj = {
        "uniprot_canonical_ac":doc["uniprot_canonical_ac"],
        "protein_name":protein_name,
        "start_pos":doc["start_pos"],
        "end_pos":doc["end_pos"],
        "site_seq":doc["site_seq"],
        "up_seq":doc["up_seq"],
        "down_seq":doc["down_seq"]
    }
    sec_list = config_obj["sitesections"]
    missing_score = 10000 
    seen = {"residue":{}, "glycosylation_type":{}, "snv_type":{}}
    for sec in sec_list:
        if sec in doc:
            record_obj[sec] = "yes" if len(doc[sec]) > 0 else "no"
            record_obj[sec+"_count"] = len(doc[sec])        
            for obj in doc[sec]:
                if "residue" in obj:
                    seen["residue"][obj["residue"]] = True
                if sec == "glycosylation":
                    gtype = obj["type"].lower() if "type" in obj else ""
                    gsubtype = obj["subtype"].lower() if "subtype" in obj else ""
                    value = "%s|%s" % (gtype, gsubtype) if gtype != "" else ""
                    seen["glycosylation_type"][value] = True
                    glytoucan_ac = obj["glytoucan_ac"]
                    if glytoucan_ac in missing_score_dict:
                        if missing_score_dict[glytoucan_ac] < missing_score:
                            missing_score = missing_score_dict[glytoucan_ac]


                if sec == "snv":
                    for value in obj["keywords"]:
                        if value.lower() in ["germline", "somatic"]:
                            seen["snv_type"][value.lower()] = True

        else:
            record_obj[sec] = ""
  
    record_obj["residue"] = seen["residue"].keys()[0] if seen["residue"] != {} else ""
    record_obj["glycosylation_type"] = ";".join(seen["glycosylation_type"].keys())
    record_obj["snv_type"] = ";".join(seen["snv_type"].keys())
    record_obj["missing_score"] = missing_score

    record_obj["tax_id"], record_obj["organism"] = "", 0
    if "species" in doc:
        if len(doc["species"]) > 0:
            record_obj["organism"] = doc["species"][0]["name"] if "name" in doc["species"][0] else ""
            record_obj["tax_id"] = doc["species"][0]["taxid"] if "taxid" in doc["species"][0] else 0


    return record_obj









def main():


    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version=" ")
    parser.add_option("-r","--recordtype",action="store",dest="recordtype",help="glycan/protein/site/motif")

    (options,args) = parser.parse_args()
    for file in ([options.recordtype]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    global config_obj
    global path_obj
    global species_obj
    global map_dict
    global data_dir
    global misc_dir
    global main_dict
    global missing_score_dict


    record_type = options.recordtype

    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
    data_dir = "reviewed/"
    misc_dir = "generated/misc/"

    file_list = glob.glob("jsondb/glycandb/*.json")
    missing_score_dict = {}
    for json_file in file_list:
        doc = json.loads(open(json_file,"r").read())
        missing_score_dict[doc["glytoucan_ac"]] = doc["missing_score"]


    #file_list = glob.glob("jsondb/%sdb/G17*.json" % (record_type))
    file_list = glob.glob("jsondb/%sdb/*P1421*.json" % (record_type))
    #file_list = glob.glob("jsondb/%sdb/Q8K4H4*.json" % (record_type))
    #file_list = glob.glob("jsondb/%sdb/GGM.000017*.json" % (record_type))
    #file_list = glob.glob("jsondb/%sdb/*.json" % (record_type))

    record_count = 0
    for json_file in file_list:
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
        if record_obj == {}:
            continue
        record_obj["record_id"] = record_id
        record_obj["record_type"] = record_type
        out_file = "jsondb/listdb/%s.json" % (record_id)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(record_obj, indent=4)))
        record_count += 1
    
    print ("make-listdb-one: final created: %s %slist records" % (record_count, record_type))





if __name__ == '__main__':
    main()

