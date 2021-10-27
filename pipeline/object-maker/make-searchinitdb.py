import os,sys
import string
import commands
from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import pymongo
from pymongo import MongoClient

sys.path.append('../../glytools/')
import libgly

__version__="1.0"
__status__ = "Dev"


def get_doc_list(coll, k_list):


    file_list = glob.glob("jsondb/" + coll.replace("c_", "") + "db/*")
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
   
    tax_id_list = []
    for species in species_list:
        tax_id_list.append(species_obj[species]["tax_id"])

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

    seen = {"glycan_type":{}, "organism":{}, "id_namespace":{}}

    k_list = [
        "glytoucan_ac", "mass", "mass_pme", "number_monosaccharides","species","composition",
        "classification", "crossref"
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


        if "number_monosaccharides" in obj:
            min_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] < min_monosaccharides else min_monosaccharides
            max_monosaccharides = obj["number_monosaccharides"] if obj["number_monosaccharides"] > max_monosaccharides else max_monosaccharides

        if "species" in obj:
            for o in obj["species"]:
                if o["taxid"] in tax_id_list:
                    org_obj_str = json.dumps({ "name": o["name"], "id":o["taxid"]})
                    seen["organism"][org_obj_str] = True
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

    res_obj["id_namespace"] = seen["id_namespace"].keys()

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
    record_type_list = res_obj.keys()

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

            d_list = tmp_dict.keys()
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
            o = {"source":d_one, "targetlist":seen[d_one].keys(), 
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
    
    doc_list = get_doc_list("c_protein", ["species", "keywords"])
    for svc in query_dict:
        for doc in doc_list:
            if query_dict[svc] not in doc["keywords"]:
                continue
            tax_id, tax_name = doc["species"][0]["taxid"], doc["species"][0]["name"]        
            if tax_id not in res_obj[svc]["organism"]:
                res_obj[svc]["organism"][tax_id] = {"id":tax_id, "name":tax_name}
            if svc == "species_to_glycoproteins":
                if "evidence_type" not in res_obj[svc]["organism"][tax_id]:
                    res_obj[svc]["organism"][tax_id]["evidence_type"] = []
                seen = {}
                for k in ["reported", "predicted"]:
                    for kw in doc["keywords"]:
                        if kw.find(k) != -1:
                            seen[k] = True
                for evdn_type in seen.keys():
                    if evdn_type not in res_obj[svc]["organism"][tax_id]["evidence_type"]:
                        res_obj[svc]["organism"][tax_id]["evidence_type"].append(evdn_type)

    svc = "species_to_glycoproteins"
    for tax_id in res_obj[svc]["organism"]:
        evdn_list = res_obj[svc]["organism"][tax_id]["evidence_type"]
        if "reported" in evdn_list and "predicted" in evdn_list:
            res_obj[svc]["organism"][tax_id]["evidence_type"].append("both")



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
        res_obj[src_record_type] = len(seen_one[src_record_type].keys())


    return res_obj



def get_protein_search_init(collection):

    tax_id_list = []
    for species in species_list: 
        tax_id_list.append(species_obj[species]["tax_id"]) 


    res_obj = {"protein_mass":{}, "organism":[], 
            "glycosylation_types":[], "id_namespace":[]}
    min_mass = 10000000.0
    max_mass = -min_mass

    seen = {"protein_mass":{}, "organism":{}, "amino_acid":{}, 
            "glycosylation_types":{}, "id_namespace":{}
    }
    k_list = ["uniprot_canonical_ac", "mass","species","glycosylation", "crossref"]

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
                        org_obj_str = json.dumps({ "name": o["name"], "id":o["taxid"]})
                        seen["organism"][org_obj_str] = True
            if "glycosylation" in obj:
                for o in obj["glycosylation"]:
                    if "residue" in o:
                        seen["amino_acid"][o["residue"]] = True
                    if "type" in o:
                        g_type = o["type"][0].upper() + o["type"][1:].lower()
                        seen["glycosylation_types"][g_type] = True
            if "crossref" in obj:
                for o in obj["crossref"]:
                    xref_badge, xref_id = o["database"], o["id"]
                    seen["id_namespace"][xref_badge] = True


    res_obj["id_namespace"] = seen["id_namespace"].keys()
    res_obj["protein_mass"]["min"] = min_mass
    res_obj["protein_mass"]["max"] = max_mass
    for org_obj_str in seen["organism"]:
        res_obj["organism"].append(json.loads(org_obj_str))
    for g_type in seen["glycosylation_types"]:
        res_obj["glycosylation_types"].append(g_type)


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

             
        config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
        aa_dict = get_aa_dict()

        species_obj = {}
        in_file = path_obj["misc"]+ "/species_info.csv"
        libgly.load_species_info(species_obj, in_file)
        
        species_list = []
        for k in species_obj:
            obj = species_obj[k]
            if obj["short_name"] not in species_list and obj["is_reference"] == "yes":
                species_list.append(obj["short_name"])

        try:

            out_obj = {}

            print "Started get_usecases_search_init ..."
            out_obj["usecases"] = get_usecases_search_init()
            print "Finished."

            print "Started get_glycan_search_init ... "
            out_obj["glycan"] = get_glycan_search_init("c_glycan")
            print "Finished."

            print "Started get_idmapping_search_init ..."
            out_obj["idmapping"] = get_idmapping_search_init()
            print "Finished."

            print "Started get_site_search_init ..."
            out_obj["site"] = get_site_search_init("c_site")
            print "Finished."

            print "Started get_super_search_init ..."
            out_obj["supersearch"] = get_super_search_init()
            print "Finished."

            print "Started get_field_limits ..."
            out_obj["fieldlimits"] = get_field_limits(config_obj["fieldlimitinfo"])
            print "Finished."

            print "Started get_protein_search_init ..."
            out_obj["protein"] = get_protein_search_init("c_protein")
            print "Finished."

            out_file = path_obj["jsondbpath"] + "/searchinitdb/searchinitdb.json"
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(out_obj, indent=4)))
            print "make-searchinitdb: final created: %s searchinitdb objects" % (1)
        except pymongo.errors.ServerSelectionTimeoutError as err:
            print err
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            print err
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}




if __name__ == '__main__':
	main()

