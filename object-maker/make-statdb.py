import os,sys
import string
from optparse import OptionParser
import glob
import json
import csvutil
import libgly

__version__="1.0"
__status__ = "Dev"



def powerset(s):
    x = len(s)
    masks = [1 << i for i in range(x)]
    for i in range(1 << x):
        yield [ss for mask, ss in zip(masks, s) if i & mask]


def get_empty_stat_obj(protein_type2kw, site_type2kw):

    # type-1 with key pattern {$param}
    # type-2 with key pattern {$param}_{$type} and $type is a value from $list
    # type-3 with key pattern {$param}_{$index_combo} and $index_combo is combo index of $list
    # type-4 with key pattern {$index_combo} and $index_combo is combo index of $list

    obj = {
        "protein":{
            "taxidlist":[]
            ,"byproteintype":{ "typelist":[]}    #type-3 with $param=tax_id, $list=typelist
            ,"bysequencetype":{"typelist":[]}   #type-2 with $param=tax_id, $list=typelist
            ,"bysitetype":{ "typelist":[]}      #type-2 with $param=tax_id, $list=typelist
            ,"glycohydrolases":{}               #type-1 with $param=tax_id
            ,"glycosyltransferases":{}          #type-1 with $param=tax_id
            ,"total":{}                         #type-1 with $param=tax_id
        }
        ,"glycan":{                             #type-4 with $list=taxidlist
            "taxidlist":[]
            ,"byglycantype":{}                  #type-1 with $param={type}_{subtype}
            ,"bymotiftype":{}                   #type-1 with $param=motif_name
            ,"bar_mass_ranges":[]
            ,"bar_sugar_ranges":[]
            ,"total":{}                         #type-1 with $param=tax_id
        }
    }

    protein_type2kw = {
        "proteins":"protein",
        "enzymes":"enzyme",
        "glycoproteins":"glycoprotein"
    }
    site_type2kw  = {
        "rwgs":"glycoprotein_reported_with_glycan",
        "rwogs":"glycoprotein_reported",
        "predicted":"glycoprotein_predicted"
    }
    
    obj["protein"]["byproteintype"]["typelist"] = list(protein_type2kw.keys())
    obj["protein"]["bysitetype"]["typelist"] = list(site_type2kw.keys())
    obj["protein"]["bysequencetype"]["typelist"] = ["canonical","isoform"]

    return obj

def get_kw_index_list(kw_list, type_list, type2kw):
    
    idx_list = []
    for t in type_list:
        kw = type2kw[t]
        if kw in kw_list:
            idx = str(type_list.index(t))
            if idx not in idx_list:
                idx_list.append(idx)
    
    #print ("kw_list:", kw_list)
    #print ("stat_fields:", type_list)
    #print ("type2kw:", type2kw)
    #print ("idx_list:", idx_list)
    #print ("\n\n")

    return idx_list


def load_species_dict(in_file):

    species2taxid = {}
    
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ref_status = row[f_list.index("is_reference")]
        if ref_status != "yes":
            continue
        tax_id = int(row[f_list.index("tax_id")])
        tax_name = row[f_list.index("glygen_name")]
        if tax_name not in species2taxid:
            species2taxid[tax_name] = tax_id

    return species2taxid


def get_old_stat_obj():


    stat_obj = {}
    taxid2species = {} 
    for species in species2taxid:
        tax_id = str(species2taxid[species])
        taxid2species[tax_id] = species
    species_map = json.loads(open("generated/misc/species_map.json", "r").read())
    tax_id_map = {}
    for tax_id in species_map:
        tax_id_map[tax_id] = species_map[tax_id]["ref_tax_id"]


    for in_file in glycan_file_list:
        doc = json.loads(open(in_file, "r").read())
        glytoucan_ac = doc["glytoucan_ac"]
        for obj in doc["species"]:
            tax_id = obj["taxid"]
            tax_id_str = str(tax_id)
            tax_id_str = tax_id_map[tax_id_str] if tax_id_str in tax_id_map else tax_id_str
            if tax_id_str not in taxid2species:
                continue
            tax_name = taxid2species[tax_id_str]
            if tax_id_str not in stat_obj:
                stat_obj[tax_id_str] = {
                    "species":tax_name,
                    "seen_glycan":{}, 
                    "seen_protein":{},
                    "seen_glycoprotein":{}
                }
            stat_obj[tax_id_str]["seen_glycan"][glytoucan_ac] = True



    for in_file in protein_file_list:
        doc = json.loads(open(in_file, "r").read())
        canon = doc["uniprot_canonical_ac"]
        for obj in doc["species"]:
            tax_id = obj["taxid"]
            tax_id_str = str(tax_id)
            tax_id_str = tax_id_map[tax_id_str] if tax_id_str in tax_id_map else tax_id_str
            if tax_id_str not in taxid2species:
                continue
            tax_name = taxid2species[tax_id_str]
            if tax_id_str not in stat_obj:
                stat_obj[tax_id_str] = {
                    "species":tax_name,
                    "seen_glycan":{},
                    "seen_protein":{},
                    "seen_glycoprotein":{}
                }
            stat_obj[tax_id_str]["seen_protein"][canon] = True
            #Glycan stats should come from glycan files only
            #for o in doc["glycosylation"]:
            #    glytoucan_ac = o["glytoucan_ac"]
            #    if glytoucan_ac != "":
            #        stat_obj[tax_id_str]["seen_glycan"][glytoucan_ac] = True
            if doc["glycosylation"] != []:
                stat_obj[tax_id_str]["seen_glycoprotein"][canon] = True


    for tax_id_str in stat_obj:
        n1 = len(list(stat_obj[tax_id_str]["seen_glycan"].keys()))
        n2 = len(list(stat_obj[tax_id_str]["seen_protein"].keys()))
        n3 = len(list(stat_obj[tax_id_str]["seen_glycoprotein"].keys()))
        stat_obj[tax_id_str]["glycans"] = n1
        stat_obj[tax_id_str]["proteins"] = n2
        stat_obj[tax_id_str]["glycoproteins"] = n3
        stat_obj[tax_id_str].pop("seen_glycan")
        stat_obj[tax_id_str].pop("seen_protein")
        stat_obj[tax_id_str].pop("seen_glycoprotein")


    return stat_obj




def get_new_stat_obj():


    protein_type2kw = { 
        "proteins":"protein",
        "enzymes":"enzyme",
        "glycoproteins":"glycoprotein"
    }
    site_type2kw  = {
        "rwgs":"glycoprotein_reported_with_glycan",
        "rwogs":"glycoprotein_reported",
        "predicted":"glycoprotein_predicted"
    }

    stat_obj = get_empty_stat_obj(protein_type2kw, site_type2kw)


    for in_file in protein_file_list:
        doc = json.loads(open(in_file, "r").read())
        for o in doc["species"]:
            tax_id, tax_name = o["taxid"], o["name"]
            tax_id_str = str(tax_id)
            if tax_id_str not in stat_obj["protein"]["taxidlist"]:
                stat_obj["protein"]["taxidlist"].append(tax_id_str)
            if tax_id_str not in stat_obj["protein"]["total"]:
                stat_obj["protein"]["total"][tax_id_str] = 0
            stat_obj["protein"]["total"][tax_id_str] += 1
            

            #--> glycosyltransferases
            if "glycosyltransferase-activity" in doc["keywords"]:
                if tax_id_str not in stat_obj["protein"]["glycosyltransferases"]:
                    stat_obj["protein"]["glycosyltransferases"][tax_id_str] = 0
                stat_obj["protein"]["glycosyltransferases"][tax_id_str] += 1
            
            #--> glycohydrolases
            if "glycohydrolase-activity" in doc["keywords"]:
                if tax_id_str not in stat_obj["protein"]["glycohydrolases"]:
                    stat_obj["protein"]["glycohydrolases"][tax_id_str] = 0
                stat_obj["protein"]["glycohydrolases"][tax_id_str] += 1
                
            #--> bysequencetype
            idx = "%s_%s" % (tax_id, "canonical")
            if idx not in stat_obj["protein"]["bysequencetype"]:
                stat_obj["protein"]["bysequencetype"][idx] = 0
            stat_obj["protein"]["bysequencetype"][idx] += 1

            idx = "%s_%s" % (tax_id, "isoform")
            if idx not in stat_obj["protein"]["bysequencetype"]:
                stat_obj["protein"]["bysequencetype"][idx] = 0
            stat_obj["protein"]["bysequencetype"][idx] += len(doc["isoforms"])


            #--> byproteintype
            type_list = stat_obj["protein"]["byproteintype"]["typelist"]
            idx_list = get_kw_index_list(doc["keywords"], type_list, protein_type2kw)
            #all of them are proteins
            if str(type_list.index("proteins")) not in idx_list: 
                idx_list.append(str(type_list.index("proteins")))
            subset_list = list(powerset(sorted(idx_list)))
            for subset in subset_list:
                if subset == []:
                    continue
                idx = "%s_%s" % (tax_id, "|".join(subset))
                if idx not in stat_obj["protein"]["byproteintype"]:
                    stat_obj["protein"]["byproteintype"][idx] = 0
                stat_obj["protein"]["byproteintype"][idx] += 1


            #--> bysitetype
            type_list = stat_obj["protein"]["bysitetype"]["typelist"]
            idx_list = get_kw_index_list(doc["keywords"], type_list, site_type2kw)
            subset_list = list(powerset(sorted(idx_list)))
            for subset in subset_list:
                if subset == []:
                    continue
                idx = "%s_%s" % (tax_id, "|".join(subset))
                if idx not in stat_obj["protein"]["bysitetype"]:
                    stat_obj["protein"]["bysitetype"][idx] = 0
                stat_obj["protein"]["bysitetype"][idx] += 1            


    mass_upper_limit_list = []
    for i in range(1, 10):
        mass_upper_limit_list.append(i*1000)
        stat_obj["glycan"]["bar_mass_ranges"].append(0)
    mass_upper_limit_list.append(1000000)
    stat_obj["glycan"]["bar_mass_ranges"].append(0)

    sugar_upper_limit_list = []
    for i in range(1, 21):
        sugar_upper_limit_list.append(i)
        stat_obj["glycan"]["bar_sugar_ranges"].append(0)
    sugar_upper_limit_list.append(1000)
    stat_obj["glycan"]["bar_sugar_ranges"].append(0)


    for in_file in glycan_file_list:
        doc = json.loads(open(in_file, "r").read())
        for i in range(0, len(sugar_upper_limit_list)):
            if "number_monosaccharides" not in doc:
                continue
            if doc["number_monosaccharides"] < sugar_upper_limit_list[i]:
                stat_obj["glycan"]["bar_sugar_ranges"][i] += 1
                break

        for i in range(0, len(mass_upper_limit_list)):
            if "mass" not in doc:
                continue
            if doc["mass"] < mass_upper_limit_list[i]:
                stat_obj["glycan"]["bar_mass_ranges"][i] += 1
                break

        # glycan by motif
        for o in doc["motifs"]:
            idx = o["name"].replace(".", ",")
            if idx not in stat_obj["glycan"]["bymotiftype"]:
                stat_obj["glycan"]["bymotiftype"][idx] = 0
            stat_obj["glycan"]["bymotiftype"][idx] += 1


        # glycan by classification
        for o in doc["classification"]:
            g_type = o["type"]["name"]
            g_subtype = o["subtype"]["name"]
            idx = "%s_%s" % (g_type, g_subtype) 
            if idx not in stat_obj["glycan"]["byglycantype"]:
                stat_obj["glycan"]["byglycantype"][idx] = 0
            stat_obj["glycan"]["byglycantype"][idx] += 1

        # glycan by species
        idx_list = []
        for o in doc["species"]:
            tax_id, tax_name = o["taxid"], o["name"]
            tax_id_str = str(tax_id)
            if tax_id_str not in stat_obj["glycan"]["taxidlist"]:
                stat_obj["glycan"]["taxidlist"].append(tax_id_str)
            idx_list.append(str(stat_obj["glycan"]["taxidlist"].index(tax_id_str)))
        if idx_list != []:
            subset_list = list(powerset(sorted(idx_list)))
            for subset in subset_list:
                if subset == []:
                    continue
                idx = "|".join(subset)
                if idx not in stat_obj["glycan"]:
                    stat_obj["glycan"][idx] = 0
                stat_obj["glycan"][idx] += 1
                if len(subset) == 1:
                    tax_id = stat_obj["glycan"]["taxidlist"][int(idx)]
                    tax_id_str = str(tax_id) 
                    if tax_id_str not in stat_obj["glycan"]["total"]:
                        stat_obj["glycan"]["total"][tax_id_str] = 0
                    stat_obj["glycan"]["total"][tax_id_str] += 1
        


    return stat_obj






###############################
def main():

    global wrk_dir
    global species2taxid
    global glycan_file_list
    global protein_file_list

    
    generated_dir = "generated/"
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    jsondb_dir = wrk_dir + "/jsondb/"


    species_file = "generated/misc/species_info.csv"
    species2taxid = load_species_dict(species_file)


    glycan_file_list = glob.glob(jsondb_dir + "/glycandb/*.json")
    protein_file_list = glob.glob(jsondb_dir + "/proteindb/*.json")

    #glycan_file_list = glob.glob(jsondb_dir + "/glycandb/G17*.json")
    #protein_file_list = glob.glob(jsondb_dir + "/proteindb/P1421*.json")
    
    old_stat_obj = get_old_stat_obj()
    #new_stat_obj = get_new_stat_obj()
    #doc = { "oldstat":old_stat_obj, "newstat":new_stat_obj}
    doc = { "oldstat":old_stat_obj}

    out_file = jsondb_dir + "/statdb/stat.json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(doc, indent=4)))


    log_file = "logs/make-statdb.log"
    msg = "make-statdb: final created: 1 objects"
    csvutil.write_log_msg(log_file, msg, "w")


            

if __name__ == '__main__':
	main()

