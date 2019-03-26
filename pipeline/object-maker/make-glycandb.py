#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict


sys.path.append('../../glytools/')
import libgly






#######################
def clean_obj(obj):

    if type(obj) is dict:
        key_list = obj.keys()
        for k1 in key_list:
            if obj[k1] in["", [], {}]:
                obj.pop(k1)
            elif type(obj[k1]) in [dict, list]:
                clean_obj(obj[k1])
    elif type(obj) is list:
        for k1 in xrange(0, len(obj)):
            if obj[k1] in["", [], {}]:
                del obj[k1]
            elif type(obj[k1]) in [dict, list]:
                clean_obj(obj[k1])
    return







#####################################
def order_obj(jsonObj):

    ordrHash = {"glytoucan_ac":1, "mass":2, "iupac":3, "wurcs":4, "glycoct":5,
            "species":6, "classification":7,"glycosylation":8,
            "enzyme":9, "crossref":10}

    for k1 in jsonObj:
        ordrHash[k1] = ordrHash[k1] if k1 in ordrHash else 1000
        if type(jsonObj[k1]) is dict:
            for k2 in jsonObj[k1]:
                ordrHash[k2] = ordrHash[k2] if k2 in ordrHash else 1000
                if type(jsonObj[k1][k2]) is dict:
                    for k3 in jsonObj[k1][k2]:
                        ordrHash[k3] = ordrHash[k3] if k3 in ordrHash else 1000
                    jsonObj[k1][k2] = OrderedDict(sorted(jsonObj[k1][k2].items(),
                        key=lambda x: float(ordrHash.get(x[0]))))
                elif type(jsonObj[k1][k2]) is list:
                    for j in xrange(0, len(jsonObj[k1][k2])):
                        if type(jsonObj[k1][k2][j]) is dict:
                            for k3 in jsonObj[k1][k2][j]:
                                ordrHash[k3] = ordrHash[k3] if k3 in ordrHash else 1000
                                jsonObj[k1][k2][j] = OrderedDict(sorted(jsonObj[k1][k2][j].items(), 
                                    key=lambda x: float(ordrHash.get(x[0]))))
            jsonObj[k1] = OrderedDict(sorted(jsonObj[k1].items(),
                key=lambda x: float(ordrHash.get(x[0]))))

    return OrderedDict(sorted(jsonObj.items(), key=lambda x: float(ordrHash.get(x[0]))))



#######################################
def main():

    config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
    root_obj =  config_obj[config_obj["server"]]["rootinfo"]

    species2taxid= {"human":9606, "mouse":10090}

    data_dir = "reviewed/"
     
    #load all dictionaries
    map_dict = {}
    dict_list_obj = json.loads(open("../../conf/glycan_dictionaries.json", "r").read())
     
    for dict_name in dict_list_obj:
        map_dict[dict_name] = {}
        ind_list = dict_list_obj[dict_name]["indexlist"]
        for pattern in dict_list_obj[dict_name]["fileglob"]:
            for in_file in glob.glob(data_dir + pattern):
                sheet_obj = {}
                libgly.load_sheet(sheet_obj, in_file, ",")
                for row in sheet_obj["data"]:
                    if row[ind_list[0]][0] == "#":
                        continue
                    key = row[ind_list[0]]
                    val = row[ind_list[1]]
                    if key not in map_dict[dict_name]:
                        map_dict[dict_name][key] = []
                    map_dict[dict_name][key].append(val)

    

    #load workbook
    work_book = {}
    file_list_obj = json.loads(open("../../conf/glycan_datasets.json", "r").read())
    for species in file_list_obj:
        for sheet_name in file_list_obj[species]:
            in_file = data_dir + "/%s_glycan_%s.csv" % (species, sheet_name)
            if sheet_name not in work_book:
                work_book[sheet_name] = {}
            libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "glytoucan_ac")

        for sheet_name in ["glycosylation_sites_unicarbkb_glytoucan"]:
            in_file = data_dir + "/%s_proteoform_%s.csv" % (species, sheet_name)
            if sheet_name not in work_book:
                work_book[sheet_name] = {}
            libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "glytoucan_ac")


    data_grid = {}
    seen_row = {}
    for sheet_name in work_book:
        for main_id in work_book[sheet_name]["data"]:
            if main_id not in data_grid:
                data_grid[main_id] = {}
            if sheet_name not in data_grid[main_id]:
                data_grid[main_id][sheet_name] = []
            for row in work_book[sheet_name]["data"][main_id]:
                s = main_id + sheet_name + json.dumps(row)
                if s not in seen_row:
                    data_grid[main_id][sheet_name].append(row)
                seen_row[s] = True



    seen = {}
    xref_obj_list = [
        {"field":"glycomedb_id", "database":"GlycomeDB", "urlkey":"glycomedb"},
        {"field":"unicarbkb_id", "database":"UniCarbKB", "urlkey":"uckb_glycan"},
        {"field":"pubchem_id", "database":"PubChem Compound", "urlkey":"pubchem_compound"},
    ]
    record_count = 0 
    doc_list = []
    for main_id in data_grid:
        doc = {}
        doc["glytoucan_ac"] = main_id
        doc["species"] = []
        species_dict = {}
        sheet_name = "idmapping"
        tmp_fl = work_book[sheet_name]["fields"]
        for tmp_row in data_grid[main_id][sheet_name]:
            tax_id = tmp_row[tmp_fl.index("taxid")]
            if tax_id not in ["9606", "10090"]:
                continue
            tax_name = map_dict["taxid2name"][tax_id][0] if tax_id in map_dict["taxid2name"] else ""
            if tax_id not in species_dict:
                species_dict[tax_id] = {"taxid":int(tax_id), "name":tax_name, "evidence":[]}
            for xref_obj in xref_obj_list:
                field,database,urlkey = xref_obj["field"], xref_obj["database"], xref_obj["urlkey"]
                field = "unicarbkb_id" if field == "uckb_id" else field
                field = "pubchem_id" if field == "pbch_id" else field
                db_id = tmp_row[tmp_fl.index(field)]
                if db_id == "":
                    continue
                if db_id[0:3] == "SID":
                    database,urlkey = "PubChem Substance", "pubchem_substance"
                db_id_clean = db_id.replace("CID", "")
                db_id_clean = db_id_clean.replace("SID", "")
                url = map_dict["xrefdb2url"][urlkey][0] % (db_id_clean)
                o = {"database":database, "id":db_id, "url":url}
                if field in ["glycomedb_id", "unicarbkb_id"]:
                    if o not in species_dict[tax_id]["evidence"]:
                        species_dict[tax_id]["evidence"].append(o)

        for tax_id in species_dict:
            doc["species"].append(species_dict[tax_id])


        doc["crossref"] = []
        seen["crossref"] = {}
        sheet_name = "idmapping"
        tmp_fl = work_book[sheet_name]["fields"]
        for tmp_row in data_grid[main_id][sheet_name]:
            for xref_obj in xref_obj_list:
                field,database,urlkey = xref_obj["field"], xref_obj["database"], xref_obj["urlkey"]
                field = "unicarbkb_id" if field == "uckb_id" else field
                field = "pubchem_id" if field == "pbch_id" else field
                db_id = tmp_row[tmp_fl.index(field)]
                if db_id == "":
                    continue
                if db_id[0:3] == "SID":
                    database,urlkey = "PubChem Substance", "pubchem_substance"
                db_id_clean = db_id.replace("CID", "")
                db_id_clean = db_id_clean.replace("SID", "")
                url = map_dict["xrefdb2url"][urlkey][0] % (db_id_clean)
                o = {"database":database, "id":db_id, "url":url}
                combo_id = "%s,%s" % (database,db_id)
                if combo_id not in seen["crossref"]:
                    doc["crossref"].append(o)
                seen["crossref"][combo_id] = True



        sheet_name = "properties"
        doc["mass"], doc["number_monosaccharides"] = -1, -1
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_row = data_grid[main_id][sheet_name][0]
                tmp_fl = work_book[sheet_name]["fields"]
                if tmp_row[tmp_fl.index("glycan_mass")] != "":
                    doc["mass"] = round(float(tmp_row[tmp_fl.index("glycan_mass")]), 2)
                if tmp_row[tmp_fl.index("monosaccharides")] != "":
                    doc["number_monosaccharides"] = int(tmp_row[tmp_fl.index("monosaccharides")])



        sheet_name = "sequences"
        doc["iupac"], doc["wurcs"], doc["glycoct"] = "", "", ""
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_row = data_grid[main_id][sheet_name][0]
                tmp_fl = work_book[sheet_name]["fields"]
                doc["iupac"] = tmp_row[tmp_fl.index("iupac_extended")]
                doc["wurcs"] = tmp_row[tmp_fl.index("wurcs")]
                doc["glycoct"] = tmp_row[tmp_fl.index("glycoct")]


        sheet_name = "motif"
        doc["motifs"] = []
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                seen["motif"] = {}
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    ac = tmp_row[tmp_fl.index("glytoucan_ac_motif")]
                    name = tmp_row[tmp_fl.index("motif_name")]
                    if ac == "":
                        continue
                    o = {"id":ac, "name":name}
                    if ac not in seen["motif"]:
                        doc["motifs"].append(o)
                    seen["motif"][ac] = True

        sheet_name = "classification"
        doc["classification"] = []
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                type_combo_list = []
                for tmp_row in data_grid[main_id][sheet_name]:
                    g_type = tmp_row[tmp_fl.index("glycan_type")].strip()
                    g_subtype = tmp_row[tmp_fl.index("glycan_subtype")].strip()
                    if g_type != "" and g_subtype != "":
                        type_combo_list.append("%s|%s" % (g_type, g_subtype))
                type_combo_list = sorted(set(type_combo_list))
                if len(type_combo_list) == 1:
                    glycan_type, glycan_subtype = type_combo_list[0].split("|")
                    type_url = map_dict["xrefdb2url"]["glycan_type"][0] % (map_dict["subtypetags"][glycan_type][0])
                    o = {"type":{"name":glycan_type, "url":type_url}}
                    if glycan_subtype != "no subtype":
                        subtype_url = map_dict["xrefdb2url"]["glycan_type"][0] % (map_dict["subtypetags"][glycan_subtype][0])
                        o["subtype"] = {"name":glycan_subtype, "url":subtype_url}
                    doc["classification"].append(o)
                else:
                    for type_combo in type_combo_list:
                        glycan_type, glycan_subtype = type_combo.split("|")
                        if glycan_subtype == "no subtype":
                            continue
                        type_url = map_dict["xrefdb2url"]["glycan_type"][0] % (map_dict["subtypetags"][glycan_type][0])
                        subtype_url = map_dict["xrefdb2url"]["glycan_type"][0] % (map_dict["subtypetags"][glycan_subtype][0])
                        o = {
                            "type":{"name":glycan_type, "url":type_url}
                            ,"subtype":{"name":glycan_subtype, "url":subtype_url}
                        }
                        doc["classification"].append(o)

        
        sheet_name = "glycosylation_sites_unicarbkb_glytoucan"
        doc["glycoprotein"] = []
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                seen_evlist = {}
                glycosylation_dict = {}
                for tmp_row in data_grid[main_id][sheet_name]:
                    pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                    canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                    uckb_id = tmp_row[tmp_fl.index("unicarbkb_id")]
                    combo_id = "%s,%s" % (canon, pos)
                    if canon not in map_dict["canon2recname"]:
                        continue
                    rec_name = map_dict["canon2recname"][canon][0] if canon in map_dict["canon2recname"] else ""
                    o = {
                        "uniprot_canonical_ac":canon,
                        "protein_name":rec_name, 
                        "position":int(pos),
                        "evidence":[]
                    }
                    if combo_id not in glycosylation_dict:
                        glycosylation_dict[combo_id] = o
                        seen_evlist[combo_id] = []
                    if uckb_id != "":
                        url = map_dict["xrefdb2url"]["uckb_glycan"][0] % (uckb_id)
                        ev_obj = {"database":"UniCarbKB", "id":uckb_id, "url":url}
                        if uckb_id not in seen_evlist[combo_id]:
                            glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                            seen_evlist[combo_id].append(uckb_id) 
                for combo_id in glycosylation_dict:
                    doc["glycoprotein"].append(glycosylation_dict[combo_id])

        sheet_name = "enzyme"
        doc["enzyme"] = []
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                seen["canon"] = {}
                for tmp_row in data_grid[main_id][sheet_name]:
                    gene_name = tmp_row[tmp_fl.index("gene_name_enzyme")]
                    canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac_enzyme")]
                    desc = map_dict["canon2recname"][canon][0] if canon in map_dict["canon2recname"] else ""
                    db_id = map_dict["canon2hgncid"][canon][0] 
                    enzyme_tax_id = tmp_row[tmp_fl.index("taxid_enzyme")]
                    gene_url = ""
                    if enzyme_tax_id == "9606":
                        gene_url = map_dict["xrefdb2url"]["hgnc"][0] % (db_id)
                    elif enzyme_tax_id == "10090":
                        gene_url = map_dict["xrefdb2url"]["mgi"][0] % (db_id)
                    o = {
                        "protein_name":desc,
                        "uniprot_canonical_ac":canon,
                        "gene":gene_name,
                        "gene_link":gene_url
                    }
                    if canon not in seen["canon"]:
                        doc["enzyme"].append(o)
                    seen["canon"][canon] = True


        doc_list.append(doc)
        record_count += 1
        if record_count%1000 == 0:
            print " ... compiled %s objects" % (record_count)
    
    print " ... final compiled %s objects" % (record_count)


    #Check image files
    img_dir = "/data/projects/glygen/downloads/glytoucan/current/images/cfg/extended/"
    for obj in doc_list:
        img_file = img_dir + "%s.png" % (obj["glytoucan_ac"])
        if os.path.isfile(img_file) == False:
            print "Warning ... file doesn't exisit:", img_file


    record_count = 0
    fout_obj = {}
    for obj in doc_list:
        cond_list= []
        cond_list.append(obj["mass"] > 0)
        cond_list.append(obj["number_monosaccharides"] > 0)
        if False not in cond_list:
            #remove empty value properties
            #clean_obj(obj)
            fout_obj[obj["glytoucan_ac"]] = order_obj(obj)
            record_count += 1
            if record_count%1000 == 0:
                print " ... filtered %s objects" % (record_count)

    out_file = path_obj["jsondbpath"] + "glycandb.json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(fout_obj, indent=4)))

    print " ... final filtered in: %s objects" % (record_count)



if __name__ == '__main__':
        main()



