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



##################
def get_sort_key_value(obj):
    return obj["date"]





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

    species_obj = config_obj["speciesinfo"]


    data_dir = "reviewed/"
    misc_dir = "generated/misc/"


    #load all dictionaries
    map_dict = {}
    dict_list_obj = json.loads(open("../../conf/glycan_dictionaries.json", "r").read())
    for dict_name in dict_list_obj:
        map_dict[dict_name] = {}
        ind_list = dict_list_obj[dict_name]["indexlist"]
        for pattern in dict_list_obj[dict_name]["fileglob"]:
            file_list = glob.glob(misc_dir + pattern)
            if len(file_list) == 0:
                print "file %s does not exist!" % (pattern)
                sys.exit()

            for in_file in glob.glob(misc_dir + pattern):
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

    #Load canon2hgncid and canon2hgncname dictionaries    
    map_dict["canon2hgncid"] = {}
    map_dict["canon2hgncname"] = {}
    file_list = glob.glob(data_dir + "/*_protein_xref_hgnc.csv") + glob.glob(data_dir + "/*_protein_xref_mgi.csv")
    for in_file in file_list:
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        for row in sheet_obj["data"]:
            key, val_one, val_two = row[0], row[1], row[2]
            if key not in map_dict["canon2hgncid"]:
                map_dict["canon2hgncid"][key] = []
            map_dict["canon2hgncid"][key].append(val_one)

            if key not in map_dict["canon2hgncname"]:
                map_dict["canon2hgncname"][key] = []
            map_dict["canon2hgncname"][key].append(val_two)
   
    #Load canon2recname dictionary
    map_dict["canon2recname"] = {}
    file_list = glob.glob(data_dir + "/*_protein_recnames.csv")
    for in_file in file_list:
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        for row in sheet_obj["data"]:
            key, val_one = row[0], row[1]
            if key not in map_dict["canon2recname"]:
                map_dict["canon2recname"][key] = []
            map_dict["canon2recname"][key].append(val_one)
    




    #load workbook
    work_book = {}
    xref_list = []
    file_list_obj = json.loads(open("../../conf/glycan_datasets.json", "r").read())
    canon2taxid = {}
    for sheet_name in file_list_obj["common"]:
        in_file = data_dir + "/glycan_%s.csv" % (sheet_name)
        if sheet_name not in work_book:
            work_book[sheet_name] = {}
        libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "glytoucan_ac")
        if sheet_name[0:5] == "xref_" and sheet_name not in xref_list:
            xref_list.append(sheet_name)
        print "Loaded %s" % (sheet_name)

    for sheet_name in ["glycosylation_sites_unicarbkb", "glycosylation_sites_harvard"]:
        if sheet_name not in work_book:
            work_book[sheet_name] = {}
        file_list = glob.glob(data_dir + "/*_proteoform_%s.csv" % (sheet_name))
        for in_file in file_list:
            libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "saccharide")
            print "Loaded %s" % (in_file)


    #make canon2taxid
    file_list = glob.glob(data_dir + "/*_protein_masterlist.csv")
    for in_file in file_list:
        species = in_file.split("/")[-1].split("_")[0]
        sheet_obj = {}
        libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "uniprotkb_canonical_ac")
        tmp_fl = sheet_obj["fields"]
        for canon in sheet_obj["data"]:
            canon2taxid[canon] = species_obj[species]["taxid"]


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
        sheet_name = "masterlist"
        doc["mass"], doc["mass_pme"], doc["number_monosaccharides"] = -1, -1, -1
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_row = data_grid[main_id][sheet_name][0]
                tmp_fl = work_book[sheet_name]["fields"]
                if tmp_row[tmp_fl.index("glycan_mass")] != "":
                    doc["mass"] = round(float(tmp_row[tmp_fl.index("glycan_mass")]), 2)
                if tmp_row[tmp_fl.index("glycan_permass")] != "":
                    doc["mass_pme"] = round(float(tmp_row[tmp_fl.index("glycan_permass")]), 2)
                if tmp_row[tmp_fl.index("monosaccharides")] != "":
                    doc["number_monosaccharides"] = int(tmp_row[tmp_fl.index("monosaccharides")])
            else:
                #if main_id not in properties, skip
                continue


        species_dict = {}
        sheet_name = "taxonomy"
        tmp_fl = work_book[sheet_name]["fields"]
        if sheet_name in data_grid[main_id]:
            for tmp_row in data_grid[main_id][sheet_name]:
                tax_id = tmp_row[tmp_fl.index("tax_id")]
                source = tmp_row[tmp_fl.index("source")]
                source_id = tmp_row[tmp_fl.index("source_id")]
                tax_name = map_dict["taxid2name"][tax_id][0] if tax_id in map_dict["taxid2name"] else ""
                if tax_id not in species_dict:
                    species_dict[tax_id] = {"taxid":int(tax_id), "name":tax_name, "evidence":[]}
                if source != "UniCarbKB":
                    url_key = source.lower()
                    url = map_dict["xrefdb2url"][url_key][0] % (source_id)
                    ev_obj = {"id":source_id, "database":source, "url":url}
                    species_dict[tax_id]["evidence"].append(ev_obj)
                else:
                    url_key = "uckb_glycan"
                    url = map_dict["xrefdb2url"][url_key][0] % (source_id)
                    ev_obj = {"id":source_id, "database":source, "url":url}
                    species_dict[tax_id]["evidence"].append(ev_obj)
        doc["species"] = []
        for tax_id in species_dict:
            doc["species"].append(species_dict[tax_id])
        

        doc["crossref"] = []
        seen["crossref"] = {}
        for sheet_name in xref_list:
            tmp_fl = work_book[sheet_name]["fields"]
            if sheet_name not in data_grid[main_id]:
                continue
            for tmp_row in data_grid[main_id][sheet_name]:
                db_id = tmp_row[tmp_fl.index("database_id")]
                db_label = tmp_row[tmp_fl.index("database_label")]
                url_key = sheet_name[5:]
                if url_key == "pubchem" and db_id[0:3] == "CID":
                    url_key = "pubchem_compound"
                if url_key == "pubchem" and db_id[0:3] == "SID":
                    url_key = "pubchem_substance"
                if url_key == "kegg":
                    url_key = "kegg_glycan"
                if url_key == "unicarbkb":
                    url_key = "uckb_glycan"
                    if db_id.find("comp_") != -1:
                        continue
                url = map_dict["xrefdb2url"][url_key][0] % (db_id)
                o = {"database":db_label, "id":db_id, "url":url}
                combo_id = "%s,%s" % (db_label,db_id)
                if combo_id not in seen["crossref"]:
                    doc["crossref"].append(o)
                    seen["crossref"][combo_id] = True


        #Extract sequences
        if main_id in data_grid:
            format2sheetname = {
                "iupac":"sequences_iupac_extended", 
                "wurcs":"sequences_wurcs", 
                "glycoct":"sequences_glycoct", 
                "inchi":"sequences_inchi", 
                "smiles_isomeric":"sequences_smiles_isomeric",
                "glycam":"sequences_glycam_iupac"
            }
            for f in format2sheetname:
                doc[f] = ""
                sheet_name = format2sheetname[f]
                if sheet_name in data_grid[main_id]:
                    tmp_row = data_grid[main_id][sheet_name][0]
                    tmp_fl = work_book[sheet_name]["fields"]
                    field_name = sheet_name.replace("sequences_", "sequence_")
                    doc[f] = tmp_row[tmp_fl.index(field_name)]

        #Extract inchi key
        doc["inchi_key"] = {}
        sheet_name = "sequences_inchi"
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    inchi_key = tmp_row[tmp_fl.index("inchi_key")]
                    url = map_dict["xrefdb2url"]["inchi_key"][0] % (inchi_key)
                    doc["inchi_key"] = {"key":inchi_key, "url":url}
        


        #Extract composition
        sheet_name = "monosaccharide_composition"
        doc["composition"] = []
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    for residue in map_dict["residue2name"]:
                        name = map_dict["residue2name"][residue][0]
                        cid = map_dict["residue2cid"][residue][0]
                        n = int(tmp_row[tmp_fl.index(residue)])
                        residue = "other" if residue.lower() == "xxx" else residue.lower()
                        url = map_dict["xrefdb2url"]["monosaccharide_residue_name"][0] % (cid)
                        o = {"name":name, "residue":residue, "count":n}
                        if cid != "":
                            o = {"name":name, "residue":residue, "count":n, "cid":cid, "url":url}
                        #if n > 0:
                        doc["composition"].append(o)

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
                    g_type = "Other" if g_type == "" else g_type
                    g_subtype = "Other" if g_subtype == "" else g_subtype
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

        #If no classification, assign other 
        if doc["classification"] == []:
            o = {"type":{"name":"Other", "url":""}, "subtype":{"name":"Other", "url":""}}
            doc["classification"].append(o)

        doc["glycoprotein"] = []
        for sheet_name in ["glycosylation_sites_unicarbkb", "glycosylation_sites_harvard"]:
            if main_id in data_grid:
                seen_evlist = {}
                glycosylation_dict = {}
                if sheet_name in data_grid[main_id]:
                    tmp_fl = work_book[sheet_name]["fields"]
                    for tmp_row in data_grid[main_id][sheet_name]:
                        pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                        canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                        amino_acid = tmp_row[tmp_fl.index("amino_acid")]
                        canon_ac = canon.split("-")[0]
                        gene_name = map_dict["canon2hgncname"][canon][0] if canon in map_dict["canon2hgncname"] else ""
                        combo_id = "%s,%s" % (canon, pos)
                        if canon not in map_dict["canon2recname"]:
                            continue
                        rec_name = map_dict["canon2recname"][canon][0] if canon in map_dict["canon2recname"] else ""
                        o = {
                            "uniprot_canonical_ac":canon,
                            "gene_name":gene_name,
                            "protein_name":rec_name, 
                            "tax_id":canon2taxid[canon],
                            "tax_name":species_obj[str(canon2taxid[canon])]["taxname"],
                            "position":int(pos),
                            "residue":amino_acid,
                            "evidence":[]
                        }
                        if combo_id not in glycosylation_dict:
                            glycosylation_dict[combo_id] = o
                            seen_evlist[combo_id] = []
                        if canon_ac != "":
                            url, src_db_id, src_db_name = "xxx", canon_ac, "xxxx"
                            if sheet_name == "glycosylation_sites_unicarbkb":
                                src_db_id = tmp_row[tmp_fl.index("unicarbkb_id")]
                                url = map_dict["xrefdb2url"]["uckb_protein"][0] % (canon_ac)
                                src_db_name = "UniCarbKB"
                            elif sheet_name == "glycosylation_sites_harvard":
                                src_db_id = tmp_row[tmp_fl.index("evidence")]
                                urlkey = "glygen_ds" if src_db_id.find("GLYDS") != -1 else "pubmed"
                                src_db_name = "GlyGen" if src_db_id.find("GLYDS") != -1 else "PubMed"
                                url = map_dict["xrefdb2url"][urlkey][0] % (src_db_id)
                            ev_obj = {"database":src_db_name, "id":src_db_id, "url":url}
                            if src_db_id not in seen_evlist[combo_id]:
                                glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                                seen_evlist[combo_id].append(src_db_id) 
                for combo_id in glycosylation_dict:
                    doc["glycoprotein"].append(glycosylation_dict[combo_id])

        sheet_name = "enzyme"
        doc["enzyme"] = []
        doc["residues"] = []
        if main_id in data_grid:
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                seen["canon"] = {}
                seen["residue"] = {}
                for tmp_row in data_grid[main_id][sheet_name]:
                    gene_name = tmp_row[tmp_fl.index("gene_name")]
                    canon = tmp_row[tmp_fl.index("uniprotkb_canonical_ac")]
                    residue_id = tmp_row[tmp_fl.index("residue_id")]
                    residue_name = tmp_row[tmp_fl.index("residue_name")]
                    parent_residue_id = tmp_row[tmp_fl.index("parent_residue_id")]

                    enzyme_tax_id = canon2taxid[canon]
                    desc = map_dict["canon2recname"][canon][0] if canon in map_dict["canon2recname"] else ""
                    db_id = map_dict["canon2hgncid"][canon][0] 
                    gene_url = ""
                    if enzyme_tax_id == "9606":
                        gene_url = map_dict["xrefdb2url"]["hgnc"][0] % (db_id)
                    elif enzyme_tax_id == "10090":
                        gene_url = map_dict["xrefdb2url"]["mgi"][0] % (db_id)

                    o = {
                        "protein_name":desc,
                        "uniprot_canonical_ac":canon,
                        "gene":gene_name,
                        "gene_link":gene_url,
                        "tax_id":canon2taxid[canon],
                        "tax_name":species_obj[str(canon2taxid[canon])]["taxname"]
                    }

                    if canon not in seen["canon"]:
                        doc["enzyme"].append(o)
                        seen["canon"][canon] = True
                    attached_by = "rxn.%s" % (canon)
                    detached_by = ""
                    r = {"id":residue_id, "name":residue_name, "attachedby":attached_by, 
                            "detachedby":detached_by, "parentid":parent_residue_id}
                    if residue_id not in seen["residue"]:
                        doc["residues"].append(r)
                        seen["residue"][residue_id] = True

        #Extract publication
        pub_dict = {}
        src_dict = {"citations_glytoucan":"GlyTouCan"}
        for sheet_name in src_dict:
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    pmid = tmp_row[tmp_fl.index("pmid")]
                    title = tmp_row[tmp_fl.index("title")]
                    journal_name = tmp_row[tmp_fl.index("journal_name")]
                    pub_date = tmp_row[tmp_fl.index("publication_date")] if "publication_date" in tmp_fl else ""
                    pub_date = pub_date.split(" ")[0]
                    source = tmp_row[tmp_fl.index("source")]
                    source_id = tmp_row[tmp_fl.index("source_id")]
                    authors = tmp_row[tmp_fl.index("authors")]
                    if pmid != "":
                        url_key = source.lower()
                        if url_key != "unicarbkb":
                            url = map_dict["xrefdb2url"][url_key][0] % (source_id)
                            ev_obj = {"id":source_id, "database":source, "url":url}
                            if pmid not in pub_dict and "" not in [pmid, title, journal_name]:
                                url = map_dict["xrefdb2url"]["pubmed"][0] % (pmid)
                                pub_obj = {"pmid":pmid, "title":title, "journal":journal_name,"url":url,
                                               "date":pub_date, "evidence":[ev_obj], "authors":authors}
                                pub_dict[pmid] = pub_obj
                            else:
                                pub_dict[pmid]["evidence"].append(ev_obj)
                        else:
                            url_key = "uckb_glycan"
                            url = map_dict["xrefdb2url"][url_key][0] % (source_id)
                            ev_obj = {"id":source_id, "database":source, "url":url}
                            if pmid not in pub_dict and "" not in [pmid, title, journal_name]:
                                url = map_dict["xrefdb2url"]["pubmed"][0] % (pmid)
                                pub_obj = {"pmid":pmid, "title":title, "journal":journal_name,"url":url,
                                            "date":pub_date, "evidence":[ev_obj], "authors":authors}
                                pub_dict[pmid] = pub_obj
                            else:
                                pub_dict[pmid]["evidence"].append(ev_obj)
            
        doc["publication"] = []
        for pmid in pub_dict:
            doc["publication"].append(pub_dict[pmid])
        doc["publication"].sort(key=get_sort_key_value, reverse=True)


        doc_list.append(doc)
        record_count += 1
        if record_count%1000 == 0:
            print " ... compiled %s objects" % (record_count)
    
    print " ... final compiled %s objects" % (record_count)


    #Check image files
    img_dir = "downloads/glytoucan/current/export/images-cfg-extended/cfg/extended/"
    for obj in doc_list:
        img_file = img_dir + "%s.png" % (obj["glytoucan_ac"])
        if os.path.isfile(img_file) == False:
            print "Warning ... no image for %s " % (obj["glytoucan_ac"])


    record_count = 0
    #fout_obj = {}
    for obj in doc_list:
        #Remove keys with -1 values
        for k in ["mass", "mass_pme", "number_monosaccharides"]:
            if obj[k] == -1:
                obj.pop(k)
        #clean_obj(obj)
        obj = order_obj(obj)
        #fout_obj[obj["glytoucan_ac"]] = order_obj(obj)

        out_file = path_obj["jsondbpath"] + "/glycandb/%s.json" % (obj["glytoucan_ac"])
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        record_count += 1
        if record_count%1000 == 0:
            print " ... created %s objects" % (record_count)
    print " ... final filtered in: %s objects" % (record_count)



    #out_file = path_obj["jsondbpath"] + "glycandb.json"
    #with open(out_file, "w") as FW:
    #    FW.write("%s\n" % (json.dumps(fout_obj, indent=4)))
    #print " ... final filtered in: %s objects" % (record_count)



if __name__ == '__main__':
        main()



