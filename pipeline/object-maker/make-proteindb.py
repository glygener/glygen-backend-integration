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

    ordrHash = {"glytoucan_ac":1, "mass":2, "iupac_extended":3, "wurcs":4, "glycoct":5,
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

    species2taxid= {"human":9606, "mouse":10090}
    taxid2name = {9606:"Homo sapiens", 10090:"Mus musculus"}

    data_dir = "reviewed/"
     
    #load all dictionaries
    map_dict = {}
    dict_list_obj = json.loads(open("../../conf/protein_dictionaries.json", "r").read())
     
    for dict_name in dict_list_obj:
        print "Loading dictionary ", dict_name
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
        


    #Load ortholog mapping
    ortholog_dict = {}
    in_file = data_dir + "/protein_ortholog.csv"
    sheet_obj = {}
    libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "uniprotkb_canonical_ac")
    tmp_fl = sheet_obj["fields"]
    for main_id in sheet_obj["data"]:
        for row in sheet_obj["data"][main_id]:
            canon_one = main_id
            canon_two = row[tmp_fl.index("uniprotkb_canonical_ac_ortholog")]
            tax_id_one = row[tmp_fl.index("tax_id")]
            tax_id_two = row[tmp_fl.index("tax_id_ortholog")]
            dataset_name = row[tmp_fl.index("evidence_dataset_name")].lower()
            dataset_id = row[tmp_fl.index("evidence_dataset_id")]
            if canon_one not in ortholog_dict:
                ortholog_dict[canon_one] = {}
            if canon_two not in ortholog_dict:
                ortholog_dict[canon_two] = {}
            if canon_two not in ortholog_dict[canon_one]:
                ortholog_dict[canon_one][canon_two] = {"taxid":tax_id_two, "evidence":[]}
            if canon_one not in ortholog_dict[canon_two]:
                ortholog_dict[canon_two][canon_one] = {"taxid":tax_id_one, "evidence":[]}
            
            url = map_dict["xrefdb2url"][dataset_name][0] % (dataset_id)
            ev_obj = {"database":dataset_name, "id":dataset_id, "url":url}
            ortholog_dict[canon_one][canon_two]["evidence"].append(ev_obj)
            ortholog_dict[canon_two][canon_one]["evidence"].append(ev_obj)




    #load workbook
    work_book = {}
    seq_hash = {}
    canon2recname = {}
    file_list_obj = json.loads(open("../../conf/protein_datasets.json", "r").read())
    #file_list_obj = json.loads(open("./tmp/protein_datasets.json", "r").read())

    for species in file_list_obj:
        for sheet_name in file_list_obj[species]:
            in_file = data_dir + "/%s_protein_%s.csv" % (species, sheet_name)
            if sheet_name not in work_book:
                work_book[sheet_name] = {}
            print "loading ", species, sheet_name
            libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "uniprotkb_canonical_ac")
        
            if sheet_name == "idmapping":
                work_book[sheet_name]["fields"].append("tax_id")
                for main_id in work_book[sheet_name]["data"]:
                    for row in work_book[sheet_name]["data"][main_id]:
                        row.append(species2taxid[species])
            if sheet_name == "recnames":
                for main_id in work_book[sheet_name]["data"]:
                    for row in work_book[sheet_name]["data"][main_id]:
                        canon2recname[main_id] = row[0]




        for sheet_name in ["glycosylation_sites_unicarbkb_glytoucan", "glycosylation_sites_uniprotkb","glycosylation_sites_pdb"]:
            in_file = data_dir + "/%s_proteoform_%s.csv" % (species, sheet_name)
            if sheet_name not in work_book:
                work_book[sheet_name] = {}
            libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "uniprotkb_canonical_ac")
    

        for sheet_name in ["allsequences"]:
            in_file = data_dir + "/%s_protein_%s.fasta" % (species, sheet_name)
            print "loading ", species, sheet_name
            for record in SeqIO.parse(in_file, "fasta"):
                seq_id = record.id.split("|")[1]
                desc = record.description
                seq_hash[seq_id] = str(record.seq.upper())


    #Load taxid mappings
    canon2taxid = {}
    for species in file_list_obj:
        in_file = data_dir + "/%s_protein_idmapping.csv" % (species)
        sheet_obj = {}                
        libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "uniprotkb_canonical_ac")
        for main_id in sheet_obj["data"]:
            canon2taxid[main_id] = species2taxid[species]

    


                        

    #Load do_id mappings
    sheet_name = "domap"
    work_book[sheet_name] = {}
    in_file = data_dir + "/human_protein_%s.csv" % (sheet_name)
    libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "do_id")
    doid2name = {}
    doid2xrefid = {}
    xrefid2doid = {}
    tmp_fl = work_book[sheet_name]["fields"]
    for do_id in work_book[sheet_name]["data"]:
        if do_id not in doid2name:
            doid2name[do_id] = []
        for tmp_row in work_book[sheet_name]["data"][do_id]:
            do_name = tmp_row[tmp_fl.index("do_name")]
            if do_name not in doid2name[do_id]:
                doid2name[do_id].append(do_name)
            
            xref_db = tmp_row[tmp_fl.index("xref_database")]
            id_in_xrefdb = tmp_row[tmp_fl.index("xref_database_id")]
            if do_id not in doid2xrefid:
                doid2xrefid[do_id] = {}
            if xref_db not in doid2xrefid[do_id]:
                doid2xrefid[do_id][xref_db] = []
            if id_in_xrefdb not in doid2xrefid[do_id][xref_db]:
                doid2xrefid[do_id][xref_db].append(id_in_xrefdb)
    
            if xref_db not in xrefid2doid:
                xrefid2doid[xref_db] = {}
            if id_in_xrefdb not in xrefid2doid[xref_db]:
                xrefid2doid[xref_db][id_in_xrefdb] = []
            if do_id not in xrefid2doid[xref_db][id_in_xrefdb]:
                xrefid2doid[xref_db][id_in_xrefdb].append(do_id)




    data_grid = {}
    seen_row = {}
    for sheet_name in work_book:
        print "transforming", sheet_name
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



    record_count = 0
    out_obj_list = []
    for main_id in data_grid:
        if "idmapping" in  data_grid[main_id]:
            #Extract uniprotkb_ac
            uniprotkb_ac = main_id.split("-")[0]

            #Extract tax_id
            tmp_fl = work_book["idmapping"]["fields"]
            tmp_row = data_grid[main_id]["idmapping"][0]


            #Extract species
            species = []            
            tax_id = canon2taxid[main_id]
            url = map_dict["xrefdb2url"]["uniprotkb"][0] % (uniprotkb_ac)
            species.append({
                "name":taxid2name[tax_id]
                ,"taxid":tax_id
                ,"evidence":[{"database":"UniProtKB", "id":uniprotkb_ac,"url":url}]
            })


            #Extract refseq accession
            refseq_ac, refseq_name, refseq_url = "", "", "" 
            sheet_name = "refseq_protein_info"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                tmp_row = data_grid[main_id][sheet_name][0]
                refseq_ac = tmp_row[tmp_fl.index("p_refseq_ac_best_match")]
                refseq_url = map_dict["xrefdb2url"]["refseq"][0] % (refseq_ac)
                refseq_name = tmp_row[tmp_fl.index("refseq_protein_name")]


            #Extract from "information"
            sheet_name = "information"
            protein_mass, protein_length, uniprotkb_id = -1.0, -1, ""
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                tmp_row = data_grid[main_id][sheet_name][0]
                protein_mass = float(tmp_row[tmp_fl.index("uniprotkb_protein_mass")])
                protein_length = int(tmp_row[tmp_fl.index("uniprotkb_protein_length")])
                uniprotkb_id = tmp_row[tmp_fl.index("uniprotkb_id")]

            #Extract from "recnames"
            sheet_name = "recnames"
            recommended_name = {}
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                tmp_row = data_grid[main_id][sheet_name][0]
                if tmp_row[tmp_fl.index("recommended_name_full")] != "":
                    recommended_name["full"] = tmp_row[tmp_fl.index("recommended_name_full")] 
                if tmp_row[tmp_fl.index("recommended_name_short")] != "":
                    recommended_name["short"] = tmp_row[tmp_fl.index("recommended_name_short")] 

            #Extract from "alternativename"
            alternative_names = []
            sheet_name = "altnames"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    full_name = tmp_row[tmp_fl.index("alternative_name_full")]
                    short_name = tmp_row[tmp_fl.index("alternative_name_short")] 
                    alternative_names += [{"full":full_name, "short":short_name}]



            #Extract "glycosylation" from unicarbkb
            glycosylation_dict = {}
            seen_glyco_ev = {}
            sheet_name = "glycosylation_sites_unicarbkb_glytoucan"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    uckb_id = tmp_row[tmp_fl.index("unicarbkb_id")]
                    pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                    glytoucan_ac =  tmp_row[tmp_fl.index("glytoucan_ac")]
                    glycosylation_type = tmp_row[tmp_fl.index("glycosylation_type")]
                    amino_acid = tmp_row[tmp_fl.index("amino_acid")].lower()
                    cond_list = []
                    cond_list.append(uckb_id != "")
                    cond_list.append(pos != "")
                    cond_list.append(amino_acid != "")
                    cond_list.append(glycosylation_type != "")
                    if False not in cond_list:
                        pos = int(pos)
                        glycosylation_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                        amino_acid = amino_acid[0].upper() + amino_acid[1:]
                        combo_id = "%s,%s,%s,%s" %(pos,glytoucan_ac,glycosylation_type,amino_acid)
                        #If this combo has not been populated
                        if combo_id not in glycosylation_dict:
                            uckb_url = map_dict["xrefdb2url"]["uckb_protein"][0] % (uniprotkb_ac)
                            glycosylation_dict[combo_id] = {
                                "glytoucan_ac":glytoucan_ac
                                ,"type":glycosylation_type
                                ,"position":pos
                                ,"residue":amino_acid
                                ,"evidence":[{"database":"UniCarbKB", "id":uckb_id, "url":uckb_url}]
                            }
                            ev_detail = "%s,%s,%s" % (combo_id,"UniCarbKB", uckb_id)
                            seen_glyco_ev[ev_detail] = True
                        pmid = tmp_row[tmp_fl.index("evidence")]
                        pmid_url = map_dict["xrefdb2url"]["pubmed"][0] % (pmid)
                        ev_obj = {"database":"PubMed", "id":pmid, "url":pmid_url}
                        ev_detail = "%s,%s,%s" % (combo_id,"PubMed", pmid)
                        if ev_detail not in seen_glyco_ev:
                            glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                            seen_glyco_ev[ev_detail] = True

            #Add glycosylation from uniprot
            sheet_name = "glycosylation_sites_uniprotkb"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                    glytoucan_ac =  ""
                    glycosylation_type = tmp_row[tmp_fl.index("glycosylation_type")]
                    amino_acid = tmp_row[tmp_fl.index("amino_acid")].lower()
                    cond_list = []
                    cond_list.append(pos != "")
                    cond_list.append(amino_acid != "")
                    cond_list.append(glycosylation_type != "")
                    if False not in cond_list:
                        pos = int(pos)
                        glytoucan_ac = ""
                        glycosylation_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                        amino_acid = amino_acid[0].upper() + amino_acid[1:]
                        combo_id = "%s,%s,%s,%s" %(pos,glytoucan_ac,glycosylation_type,amino_acid)
                        uniprot_url = map_dict["xrefdb2url"]["uniprotkb"][0] % (uniprotkb_ac + "#ptm_processing")
                        ev_obj = {"database":"UniProtKB", "id":main_id, "url":uniprot_url}
                        ev_detail = "%s,%s,%s" % (combo_id,"UniProtKB", main_id)
                        #If this combo has not been populated 
                        if combo_id not in glycosylation_dict:
                            glycosylation_dict[combo_id] = {
                                "glytoucan_ac":glytoucan_ac
                                ,"type":glycosylation_type
                                ,"position":pos
                                ,"residue":amino_acid
                                ,"evidence":[ev_obj]
                            }
                        elif ev_detail not in seen_glyco_ev:
                            #In case this combo has been populated, add UniProtKB evidence
                            glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                            seen_glyco_ev[ev_detail] = True

                        #Add pmid to evidence bucket
                        if tmp_row[tmp_fl.index("data_source")] == "PubMed":
                            pmid = tmp_row[tmp_fl.index("evidence")]
                            pmid_url = map_dict["xrefdb2url"]["pubmed"][0] % (pmid)
                            ev_obj = {"database":"PubMed", "id":pmid, "url":pmid_url}
                            ev_detail = "%s,%s,%s" % (combo_id,"PubMed", pmid)
                            if ev_detail not in seen_glyco_ev:
                                glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                                seen_glyco_ev[ev_detail] = True
            

            #Add glycosylation from pdb
            sheet_name = "glycosylation_sites_pdb"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    pdb_id =  tmp_row[tmp_fl.index("pdb_id")]
                    pos = tmp_row[tmp_fl.index("glycosylation_site_uniprotkb")]
                    pdb_pos = tmp_row[tmp_fl.index("glycosylation_site_pdb")]
                    glytoucan_ac =  tmp_row[tmp_fl.index("glytoucan_ac")]
                    glycosylation_type = tmp_row[tmp_fl.index("glycosylation_type")]
                    amino_acid = tmp_row[tmp_fl.index("amino_acid")].lower()
                    cond_list = []
                    cond_list.append(pos != "")
                    cond_list.append(amino_acid != "")
                    cond_list.append(glycosylation_type != "")
                    if False not in cond_list:
                        pos = int(pos)
                        pdb_pos = int(pdb_pos)
                        glytoucan_ac = ""
                        glycosylation_type = obj["glycosylation_type"][0]
                        glycosylation_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                        amino_acid = amino_acid[0].upper() + amino_acid[1:]
                        combo_id = "%s,%s,%s,%s" %(pos,glytoucan_ac,glycosylation_type,amino_acid)      
                        pdb_url = map_dict["xrefdb2url"]["pdb4glycosylation"][0] % (pdb_id)
                        pdb_id_lbl = pdb_id + " (position %s)" % (pdb_pos)
                        ev_obj = {"database":"PDB", "id":pdb_id_lbl, "url":pdb_url}
                        ev_detail = "%s,%s,%s" % (combo_id,"PDB", pdb_id)
                        #If this combo has not been populated
                        if combo_id not in glycosylation_dict:
                            glycosylation_dict[combo_id] = {
                                "glytoucan_ac":glytoucan_ac
                                ,"type":glycosylation_type
                                ,"position":pos
                                ,"residue":amino_acid
                                ,"evidence":[ev_obj]
                            }
                        elif ev_detail not in seen_glyco_ev:
                            #If this combo has been populated, add PDB evidence
                            glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                            seen_glyco_ev[ev_detail] = True

            #Now make glycosylation list
            glycosylation_list = []
            for combo_id in glycosylation_dict:
                glycosylation_list.append(glycosylation_dict[combo_id])


            #Extract from "pathway"
            pathway = []
            sheet_name = "xref_kegg"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    db_id = tmp_row[tmp_fl.index("database_id")]
                    db_label = tmp_row[tmp_fl.index("database_label")]
                    url = map_dict["xrefdb2url"]["kegg"][0] % (db_id)
                    pathway.append({"id":db_id, "name":db_label, "resource":"KEGG", "url":url})

            sheet_name = "xref_reactome"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    db_id = tmp_row[tmp_fl.index("database_id")]
                    db_label = tmp_row[tmp_fl.index("database_label")]
                    url = map_dict["xrefdb2url"]["reactome"][0] % (db_id)
                    pathway.append({"id":db_id, "name":db_label, "resource":"Reactome", "url":url})


            #Extract enzyme roles (glycosyltransferases, glycohydrolases)
            keywords = []
            sheet_name = "glycosyltransferase"
            if sheet_name in data_grid[main_id]:
                keywords.append("glycosyltransferase-activity")
            sheet_name = "glycohydrolase"
            if sheet_name in data_grid[main_id]:
                keywords.append("glycohydrolase-activity")


            #Extract function
            function_list = []
            sheet_name = "function_uniprot"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                seen = {"ann":{}}
                for tmp_row in data_grid[main_id][sheet_name]:
                    ann = tmp_row[tmp_fl.index("annotation")]
                    if ann not in seen["ann"]:
                        database_name = tmp_row[tmp_fl.index("database_name")]
                        database_id = tmp_row[tmp_fl.index("database_id")]
                        database_id = database_id.split("-")[0]
                        database_url = map_dict["xrefdb2url"]["uniprotkb"][0] % (database_id + "#function")
                        ev_list = [{"database":database_name, "id":database_id, "url":database_url}]
                        function_obj = { "annotation":ann, "url":database_url, "evidence":ev_list}
                        function_list.append(function_obj)
                        seen["ann"][ann] = True


            sheet_name = "function_refseq"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                seen = {"ann":{}}
                for tmp_row in data_grid[main_id][sheet_name]:
                    ann = tmp_row[tmp_fl.index("annotation")]
                    if ann not in seen["ann"]:
                        database_name = tmp_row[tmp_fl.index("database_name")]
                        database_id = tmp_row[tmp_fl.index("database_id")]
                        database_url = map_dict["xrefdb2url"]["refseq"][0] % (database_id)
                        ev_list = [{"database":database_name, "id":database_id, "url":database_url}]
                        function_obj = { "annotation":ann, "url":database_url, "evidence":ev_list}
                        function_list.append(function_obj)
                        seen["ann"][ann] = True

            #xxxxx
            #tmp_fl = work_book[sheet_name]["fields"]
            #for tmp_row in data_grid[main_id][sheet_name]:
            #    val = tmp_row[tmp_fl.index("")]



            #Extract disease
            doid2icd10 = {}
            disease_dict = {}
            sheet_name = "disease"
            seen = {"omim":{}}
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    database = tmp_row[tmp_fl.index("database_name")]
                    id_in_database = tmp_row[tmp_fl.index("database_id")]
                    label_in_database = tmp_row[tmp_fl.index("database_label")]
                    if database == "omim":
                        if id_in_database in xrefid2doid["omim"]:
                            url = map_dict["xrefdb2url"]["omim"][0] % (id_in_database)
                            ev_obj = {"database":"OMIM", "id":id_in_database,"url":url}
                            for do_id in xrefid2doid["omim"][id_in_database]:
                                do_name = doid2name[do_id][0]
                                do_name = do_name[0].upper() + do_name[1:]
                                icd10_id = ""
                                if "icd10cmid" in doid2xrefid[do_id]:
                                    icd10_id = doid2xrefid[do_id]["icd10cmid"][0]
                                url = map_dict["xrefdb2url"]["do"][0] % (do_id)
                                if do_id not in disease_dict:
                                    o = {"name":do_name,"doid":do_id,"icd10":icd10_id,"url":url, "evidence":[]}
                                    disease_dict[do_id] = o
                                if id_in_database not in seen["omim"]:
                                    disease_dict[do_id]["evidence"].append(ev_obj)
                                    seen["omim"][id_in_database] = True

            disease_list = []
            for do_id in disease_dict:
                disease_list.append(disease_dict[do_id])


            #Extract mutation
            mutation_list = []
            sheet_name = "mutation"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    start_pos = int(tmp_row[tmp_fl.index("aa_pos")])
                    end_pos = start_pos
                    ref_aa = tmp_row[tmp_fl.index("ref_aa")]
                    alt_aa = tmp_row[tmp_fl.index("alt_aa")]
                    freq = tmp_row[tmp_fl.index("mut_freq")]
                    do_id = tmp_row[tmp_fl.index("do_id")]
                    do_name = tmp_row[tmp_fl.index("do_name")].split("/")[1].strip()
                    do_name = do_name[0].upper() + do_name[1:]
                    s = "High frequency mutation in"
                    ann = "%s %s (DOID:%s) observed in %s subjects" % (s, do_name,do_id,freq)
                    icd10_id = ""
                    database_url = map_dict["xrefdb2url"]["do"][0] % (do_id)
                    disease_obj = {"name":do_name, "icd10":icd10_id, "doid":do_id, "url":database_url}
                    biomuta_url = map_dict["xrefdb2url"]["biomuta"][0] % (uniprotkb_ac)
                    ev_list = [{"database":"BioMuta", "id":uniprotkb_ac, "url":biomuta_url}]
                    mutation_obj = {
                        "annotation":ann, "type":"Point mutation", "start_pos":start_pos, "end_pos":end_pos,
                        "sequence_org":ref_aa, "sequence_mut":alt_aa,  
                        "disease":disease_obj, 
                        "evidence":ev_list
                    }
                    mutation_list.append(mutation_obj)



            #Extract orthologs
            orthologs = []
            if main_id in ortholog_dict:
                for ortholog_canon in ortholog_dict[main_id]:
                    tax_id = int(ortholog_dict[main_id][ortholog_canon]["taxid"])
                    org_name = taxid2name[tax_id]
                    recname = canon2recname[ortholog_canon] if ortholog_canon in canon2recname else ""
                    seq = seq_hash[ortholog_canon]
                    obj = {
                        "uniprot_canonical_ac":ortholog_canon
                        ,"protein_name":recname
                        ,"tax_id":tax_id
                        ,"organism":org_name
                        ,"evidence":[]
                        ,"sequence": {
                            "sequence": seq
                            ,"length": len(seq)
                        }
                    }
                    for o in ortholog_dict[main_id][ortholog_canon]["evidence"]:
                        obj["evidence"].append(o)
                    orthologs.append(obj)


            #Extract transcript
            isoform2locusobj = {}
            sheet_name = "transcriptlocus"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    isoform = tmp_row[tmp_fl.index("uniprotkb_isoform_ac")]
                    transcript_id = tmp_row[tmp_fl.index("transcript_id")]
                    peptide_id = tmp_row[tmp_fl.index("peptide_id")]
                    transcript_url = map_dict["xrefdb2url"]["ensembl"][0] % (transcript_id) 
                    peptide_url = map_dict["xrefdb2url"]["ensembl"][0] % (peptide_id)
                    start_pos = int(tmp_row[tmp_fl.index("start_pos")])
                    end_pos = int(tmp_row[tmp_fl.index("end_pos")])
                    chr_id = tmp_row[tmp_fl.index("chromosome_id")]
                    isoform2locusobj[isoform] = {
                        "chromosome":chr_id
                        ,"start_pos":start_pos
                        ,"end_pos":end_pos
                        ,"evidence":[
                            {
                                "database":"Ensembl Transcript"
                                ,"id":transcript_id
                                ,"url":transcript_url
                            }
                            ,{
                                "database":"Ensembl Peptide"
                                ,"id":peptide_id
                                 ,"url":peptide_url
                            }
                        ]
                    }
            
            #Extract transcript
            genelocus_dict = {}
            sheet_name = "genelocus"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                tmp_row = data_grid[main_id][sheet_name][0]
                ensemble_gene_id = tmp_row[tmp_fl.index("ensembl_gene_id")]
                start_pos = int(tmp_row[tmp_fl.index("start_pos")])
                end_pos = int(tmp_row[tmp_fl.index("end_pos")])
                chr_id = tmp_row[tmp_fl.index("chromosome_id")]
                gene_url = map_dict["xrefdb2url"]["ensembl"][0] % (ensemble_gene_id)
                genelocus_dict[main_id] = {
                    "chromosome":chr_id
                    ,"start_pos":start_pos
                    ,"end_pos":end_pos
                    ,"evidence":[
                        {
                            "database":"Ensembl Gene"
                            ,"id":ensemble_gene_id
                            ,"url":gene_url
                        }
                    ]
                }

            crossref = []
            #Add pdb references
            bgee_id = ""
            gene = []
            for xrefdb in ["pdb", "kegg", "reactome", "bgee", "hgnc", "mgi"]:
                sheet_name = "xref_" + xrefdb
                db_name = config_obj["xref"][sheet_name]["dbname"]
                if sheet_name in data_grid[main_id]:
                    tmp_fl = work_book[sheet_name]["fields"]
                    for tmp_row in data_grid[main_id][sheet_name]:
                        db_label = tmp_row[tmp_fl.index("database_label")]
                        db_id = tmp_row[tmp_fl.index("database_id")]
                        db_url = map_dict["xrefdb2url"][xrefdb][0] % (db_id)
                        o = {"database":db_name, "id":db_id, "url":db_url}
                        crossref.append(o)
                        if xrefdb == "bgee":
                            bgee_id = db_id
                        if xrefdb in ["hgnc", "mgi"] and db_label != "":
                            locus_obj = genelocus_dict[main_id] if main_id in genelocus_dict else {}
                            gene.append({"name":db_label, "url":db_url, "locus":locus_obj})


            #Extract expression_normal
            seen_doid = {}
            seen_uberonid = {}
            expression_tissue_list = []
            sheet_name = "expression_normal"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    uberon_id = tmp_row[tmp_fl.index("uberon_anatomy_id")].split(":")[1]
                    if uberon_id in seen_uberonid:
                        continue
                    seen_uberonid[uberon_id] = True
                    uberon_name = tmp_row[tmp_fl.index("uberon_name")]
                    uberon_name = uberon_name[0].upper() + uberon_name[1:]
                    call = "yes" if tmp_row[tmp_fl.index("expression_call")] == "present" else "no"
                    url = map_dict["xrefdb2url"]["uberon"][0] % (uberon_id)
                    tissue_obj = {"name":uberon_name,"uberon":uberon_id, "url":url}
                    url = map_dict["xrefdb2url"]["bgee"][0] % (bgee_id)
                    ev_list = [{"database":"Bgee", "id":bgee_id, "url":url}]
                    expression_tissue_list.append({"tissue":tissue_obj, "present":call, "evidence":ev_list})
                    if uberon_id in map_dict["uberonid2doid"]:
                        for do_id in map_dict["uberonid2doid"][uberon_id]:
                            seen_doid[do_id] = True


            #Extract expression_disease
            expression_disease_list = []
            sheet_name = "expression_disease"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    do_id = tmp_row[tmp_fl.index("do_id")]
                    do_name = tmp_row[tmp_fl.index("do_name")]
                    parent_doid = tmp_row[tmp_fl.index("parent_doid")]
                    significance = tmp_row[tmp_fl.index("significance")]
                    direction = tmp_row[tmp_fl.index("direction")]
                    #disregard if both doid and parent_doi are not mapped to uberonid
                    if do_id not in seen_doid and parent_doid not in seen_doid:
                        continue
                    icd10 = doid2xrefid[do_id]["icd10cmid"][0] if "icd10cmid" in doid2xrefid[do_id] else ""
                    url = map_dict["xrefdb2url"]["do"][0] % (do_id)
                    disease_obj = {"name":do_name,"doid":do_id, "icd10":icd10, "url":url}
                    url = map_dict["xrefdb2url"]["bioxpress"][0] % (uniprotkb_ac)
                    ev_list = [{"database":"BioXpress", "id":uniprotkb_ac, "url":url}]
                    expression_disease_list.append({"disease":disease_obj, 
                        "significant":significance, "trend":direction,"evidence":ev_list})
                                                                            

            #Extract sequence for canonical
            sequence = {
                "sequence":seq_hash[main_id]
                ,"length":len(seq_hash[main_id])
            }


            #Extract isoforms
            isoforms = []
            sheet_name = "idmapping"
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    reviewed_isoform = tmp_row[tmp_fl.index("reviewed_isoforms")]
                    unreviewed_isoform = tmp_row[tmp_fl.index("unreviewed_isoforms")]
                    isoform_list = []
                    isoform_list += [reviewed_isoform] if reviewed_isoform != "" else []
                    isoform_list += [unreviewed_isoform] if unreviewed_isoform != "" else []
                    for isoform in isoform_list:
                        isoform_seq = seq_hash[isoform] if isoform in seq_hash else ""
                        locus_obj = {
                            "chromosome":"",
                            "start_pos":0, 
                            "end_pos":0,
                            "evidence":[]
                        }
                        locus_obj = isoform2locusobj[isoform] if isoform in isoform2locusobj else locus_obj
                        isoforms.append({
                            "isoform_ac":isoform
                            ,"url":map_dict["xrefdb2url"]["uniprot_isoform"][0] % (uniprotkb_ac, isoform)
                            ,"sequence":{
                                "sequence":isoform_seq
                                ,"length":len(isoform_seq)
                            },
                            "locus":locus_obj
                        })


            #Extract publication
            publication = []
            sheet_name = "citations"
            seen_publication = {}
            if sheet_name in data_grid[main_id]:
                tmp_fl = work_book[sheet_name]["fields"]
                for tmp_row in data_grid[main_id][sheet_name]:
                    pmid = tmp_row[tmp_fl.index("pmid")]
                    title = tmp_row[tmp_fl.index("title")]
                    journal_name = tmp_row[tmp_fl.index("journal_name")]
                    if pmid != "":
                        url = map_dict["xrefdb2url"]["pubmed"][0] % (pmid)
                        if pmid not in seen_publication and "" not in [pmid, title, journal_name]:
                            pub_obj = {"pmid":pmid, "title":title, "journal":journal_name,"url":url}
                            publication.append(pub_obj)
                            seen_publication[pmid] = True

            out_obj = {
                    "uniprot_canonical_ac":main_id
                    ,"uniprot_id":uniprotkb_id
                    ,"refseq":{
                        "ac":refseq_ac
                        ,"name":refseq_name
                        ,"url":refseq_url

                    }
                    ,"mass":{
                        "chemical_mass":protein_mass
                        ,"monoisotopic_mass":protein_mass
                    }
                    ,"function":function_list
                    ,"disease":disease_list
                    ,"mutation":mutation_list
                    ,"expression_tissue":expression_tissue_list
                    ,"expression_disease":expression_disease_list
                    ,"recommendedname":recommended_name
                    ,"alternativenames":alternative_names
                    ,"species":species
                    ,"gene":gene
                    ,"crossref":crossref
                    ,"sequence":sequence
                    ,"isoforms":isoforms
                    ,"glycosylation":glycosylation_list
                    ,"keywords":keywords
                    ,"publication":publication
                    ,"pathway":pathway
                    ,"orthologs":orthologs
            }
            out_obj_list.append(out_obj)
            record_count += 1
            if record_count%1000 == 0:
                print " ... compiled %s objects" % (record_count)

    
    print " ... final compiled %s objects" % (record_count)


    #make dictionary of protein names to name ortholog protein 
    canon2protein_name = {}
    for obj in out_obj_list:
        uniprot_canonical_ac = obj["uniprot_canonical_ac"]
        protein_name = obj["recommendedname"]["full"] if "full" in obj["recommendedname"] else ""
        canon2protein_name[uniprot_canonical_ac] = protein_name
        
    for obj in out_obj_list:
        for o in obj["orthologs"]:
            uniprot_canonical_ac = o["uniprot_canonical_ac"]
            if "idmapping" in  data_grid[uniprot_canonical_ac]:
                o["protein_name"] = canon2protein_name[uniprot_canonical_ac]
            else:
                o["protein_name"] = ""


    fout_obj = {}
    record_count = 0
    for obj in out_obj_list:
        cond_list= []
        if False not in cond_list:
            #clean_obj(obj)
            fout_obj[obj["uniprot_canonical_ac"]] = order_obj(obj)
            record_count += 1 
            if record_count%1000 == 0:
                print " ... filtered %s objects" % (record_count)
        
    out_file = path_obj["jsondbpath"] + "proteindb.json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(fout_obj, indent=4)))
    print " ... final filtered in: %s objects" % (record_count)




if __name__ == '__main__':
        main()



