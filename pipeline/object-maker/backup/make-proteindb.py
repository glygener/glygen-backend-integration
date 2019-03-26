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
import random



__version__="1.0"
__status__ = "Dev"



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



    ordrHash = {"uniprot_canonical_ac":1, "uniprot_canonical_id":2, "mass":3

    }

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


#####################
def random_string(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))




#######################################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
	root_obj =  config_obj[config_obj["server"]]["rootinfo"]


	db_obj = {}
        species2taxid= {"human":9606, "mouse":10090}
        taxid2name = {9606:"Homo sapiens", 10090:"Mus musculus"}
		
        #data_dir = path_obj["htmlpath"] + "/csv/"
        data_dir = "/data/projects/glygen/generated/datasets/reviewed/"

        file_list_obj = json.loads(open(data_dir + "/protein_file_list.json", "r").read())
        file_list1, file_list2 = [], []
        for file_name_part in file_list_obj["protein"]:
            file_list1.append(data_dir + "/human_protein_" + file_name_part + ".csv")
            file_list1.append(data_dir + "/mouse_protein_" + file_name_part+ ".csv")
        for file_name_part in file_list_obj["proteoform"]:
            file_list2.append(data_dir + "/human_proteoform_" + file_name_part + ".csv")
            file_list2.append(data_dir + "/mouse_proteoform_" + file_name_part + ".csv")




        fasta_files = glob.glob(data_dir + "*_protein_all.fasta")
        seq_hash = {}
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = record.id.split("|")[1]
                desc = record.description
                seq_hash[seq_id] = str(record.seq.upper())
       

        uberonid2doidlist = {}
        map_file = data_dir + "/doid2uberonid-mapping.csv"
        with open(map_file, 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                row_count += 1
                if row_count == 1:
                    continue
                if row[1] not in uberonid2doidlist:
                    uberonid2doidlist[row[1]] = []
                uberonid2doidlist[row[1]].append(row[0].strip())
    
        aa_format = {"full2three":{}, "full2one":{}}
        map_file = data_dir + "/amino_acid_format.csv"
        with open(map_file, 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                row_count += 1
                if row_count == 1:
                    continue
                aa_format["full2three"][row[0].lower()] = row[1].strip().lower()
                aa_format["full2one"][row[0].lower()] = row[2].strip().lower()






        url_dict = {}
        url_templates_file = data_dir + "/database_url_templates.csv"
        with open(url_templates_file, 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                row_count += 1
                if row_count == 1:
                    continue
                if row[0][1] == "#":
                    continue
                url_dict[row[0]] = row[1].strip()


        FW = open("tmp/protein-fields.txt", "w")
        for in_file in file_list1:
                file_name = os.path.basename(in_file)
                if file_name in file_list_obj["exceptionlist"]:
                    print " ... exception! %s is ignored! " % (file_name)
                    continue
                if os.path.isfile(in_file) == False:
                    print " ... warning! %s does not exisit! " % (file_name)
                    continue
                prefix = "_".join(file_name.split(".")[0].split("_")[2:])
                org_name = file_name.split("_")[0]
                tax_id = species2taxid[org_name]
                with open(in_file, 'r') as FR:
			dataFrame = csv.reader(FR, delimiter=',', quotechar='"')
			rowCount = 0
			field_list = []
			for row in dataFrame:
                            rowCount += 1
			    if rowCount == 1:
                                field_list = row
                                FW.write(">%s %s\n%s\n\n" % (file_name,prefix,"\n".join(row)))
                            else:
                                main_id_index = -1
                                if "uniprotkb_acc_canonical" in field_list:
                                    main_id_index = field_list.index("uniprotkb_acc_canonical")
                                elif "uniprotkb_canonical_ac" in field_list:
                                    main_id_index = field_list.index("uniprotkb_canonical_ac")
                                else:
                                    main_id_index = field_list.index("uniprot_canonical_ac")
                                id1 = row[main_id_index]
				if id1 not in db_obj:
                                    db_obj[id1] = {"taxid":tax_id}
				if prefix not in db_obj[id1]:
				    db_obj[id1][prefix] = []
				row_obj = {}
                                for j in xrange(1,len(field_list)):
                                    id2 = field_list[j]
                                    row_obj[id2] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                                    if "" in row_obj[id2]:
                                        row_obj[id2].remove("")
                                db_obj[id1][prefix].append(row_obj)
                print " ...... done! %s (tax_id=%s)" % (file_name, tax_id)


        for in_file in file_list2:
                file_name = os.path.basename(in_file)
                if file_name in file_list_obj["exceptionlist"]:
                    print " ... exception! %s is ignored! " % (file_name)
                    continue
                if os.path.isfile(in_file) == False:
                    print " ... warning! %s does not exisit! " % (file_name)
                    continue
                prefix = "_".join(file_name.split(".")[0].split("_")[2:])
                org_name = file_name.split("_")[0]
                tax_id = species2taxid[org_name]
                with open(in_file, 'r') as FR:
                        dataFrame = csv.reader(FR, delimiter=',', quotechar='"')
                        rowCount = 0
                        field_list = []
                        for row in dataFrame:
                                rowCount += 1
                                if rowCount == 1:
                                        field_list = row
                                        FW.write(">%s %s\n%s\n\n" % (file_name,prefix,"\n".join(row)))
                                else:
                                        main_id_index = -1
                                        if "uniprotkb_acc_canonical" in field_list:
                                            main_id_index = field_list.index("uniprotkb_acc_canonical")
                                        else:
                                            main_id_index = field_list.index("uniprot_canonical_ac")
                                        id1 = row[main_id_index]
                                        if id1 not in db_obj:
                                                db_obj[id1] = {"taxid":tax_id}
                                        if prefix not in db_obj[id1]:
                                                db_obj[id1][prefix] = []
                                        row_obj = {}
                                        for j in xrange(0,len(field_list)):
                                                id2 = field_list[j]
                                                row_obj[id2] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                                                if "" in row_obj[id2]:
                                                    row_obj[id2].remove("") 
                                        db_obj[id1][prefix].append(row_obj)
                print " ...... done! %s (tax_id=%s)" % (file_name, tax_id)
        FW.close()

        record_count = 0
        out_obj_list = []
	for id1 in db_obj:
            if "idmapping" in  db_obj[id1]:

                #Extract accession
                uniprotkb_ac = id1.split("-")[0]
                tax_id = db_obj[id1]["taxid"]
                #Extract refseq accession
                refseq_ac, refseq_name, refseq_url = "", "", "" 
                main_key = "refseq_protein_info"
                if main_key in db_obj[id1]:
                    if len(db_obj[id1][main_key][0]["p_refseq_acc_best_match"]) > 0:
                        refseq_ac = db_obj[id1][main_key][0]["p_refseq_acc_best_match"][0]
                        refseq_url = url_dict["refseq"] % (refseq_ac)
                        refseq_name = db_obj[id1][main_key][0]["refseq_protein_name"][0]


                #Extract from "information"
                main_key = "information"
                species = []
                keywords = []
                protein_mass, protein_length, uniprotkb_id = -1.0, -1, ""
                if main_key in db_obj[id1]:
                    #if len(db_obj[id1][main_key][0]["keywords_label"]) > 0:
                    #    keywords = db_obj[id1][main_key][0]["keywords_label"]
                    if len(db_obj[id1][main_key][0]["protein_mass"]) > 0:
                        protein_mass = float(db_obj[id1][main_key][0]["protein_mass"][0])
                    if len(db_obj[id1][main_key][0]["protein_length"]) > 0:
                        protein_length = int(db_obj[id1][main_key][0]["protein_length"][0])
                    url = url_dict["uniprotkb"] % (uniprotkb_ac)
                    species.append({
                            "name":taxid2name[tax_id]
                            ,"taxid":tax_id
                            ,"evidence":[{"database":"UniProtKB", "id":uniprotkb_ac,"url":url}]
                        }
                    )
                    if len(db_obj[id1][main_key][0]["uniprotkb_id"]) > 0:
                        uniprotkb_id = db_obj[id1][main_key][0]["uniprotkb_id"][0]
                
                
                #Extract from "recnames"
                main_key = "recnames"
                recommended_name = {}
                if main_key in db_obj[id1]:
                    if len(db_obj[id1][main_key][0]["recommended_name_full"]) > 0:
                        recommended_name["full"] = db_obj[id1][main_key][0]["recommended_name_full"][0]
                    if len(db_obj[id1][main_key][0]["recommended_name_short"]) > 0:
                        recommended_name["short"] = db_obj[id1][main_key][0]["recommended_name_short"][0]

                #Extract from "alternativename"
                alternative_names = []
                main_key = "alternativename"
                if main_key in db_obj[id1]:
                    for o in db_obj[id1][main_key]:
                        full_name, short_name = "", ""
                        if len(o["alternative_name_full"]) > 0:
                            full_name = o["alternative_name_full"][0]
                        if len(o["alternative_name_short"]) > 0:
                            short_name = o["alternative_name_short"][0]
                        alternative_names += [{"full":full_name, "short":short_name}]


                
                #Extract "glycosylation" from unicarbkb
                glycosylation_dict = {}
                seen_glyco_ev = {}
                main_key = "glycosylation_sites_unicarbkb_glytoucan"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        cond_list = []
                        cond_list.append(len(obj["uckb_id"]) > 0)
                        cond_list.append(len(obj["glycosylation_site"]) > 0)
                        cond_list.append(len(obj["amino_acid"]) > 0)
                        cond_list.append(len(obj["glycosylation_type"]) > 0)
                        if False not in cond_list:
                            pos = int(obj["glycosylation_site"][0])
                            glytoucan_ac = obj["glytoucan_acc"][0] if len(obj["glytoucan_acc"]) > 0 else ""
                            glycosylation_type = obj["glycosylation_type"][0]
                            glycosylation_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                            amino_acid = obj["amino_acid"][0].lower()
                            amino_acid = amino_acid[0].upper() + amino_acid[1:]
                            combo_id = "%s,%s,%s,%s" %(pos,glytoucan_ac,glycosylation_type,amino_acid)
                            #If this combo has not been populated
                            if combo_id not in glycosylation_dict:
                                uckb_id = obj["uckb_id"][0]
                                uckb_url = url_dict["uckb_protein"] % (uniprotkb_ac)
                                glycosylation_dict[combo_id] = {
                                    "glytoucan_ac":glytoucan_ac
                                    ,"type":glycosylation_type
                                    ,"position":pos
                                    ,"residue":amino_acid
                                    ,"evidence":[{"database":"UniCarbKB", "id":uckb_id, "url":uckb_url}]
                                }
                                ev_detail = "%s,%s,%s" % (combo_id,"UniCarbKB", uckb_id)
                                seen_glyco_ev[ev_detail] = True

                            #Add pmid to evidence bucket
                            for ev in obj["evidence"]:
                                pmid = ev
                                pmid_url = url_dict["pubmed"] % (pmid)
                                ev_obj = {"database":"PubMed", "id":pmid, "url":pmid_url}
                                ev_detail = "%s,%s,%s" % (combo_id,"PubMed", pmid)
                                if ev_detail not in seen_glyco_ev:
                                    glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                                    seen_glyco_ev[ev_detail] = True


                #Add glycosylation from uniprot
                main_key = "glycosylation_sites_uniprot"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        cond_list = []
                        cond_list.append(len(obj["glycosylation_site"]) > 0)
                        cond_list.append(len(obj["amino_acid"]) > 0)
                        cond_list.append(len(obj["glycosylation_type"]) > 0)
                        if False not in cond_list:
                            pos = int(obj["glycosylation_site"][0])
                            glytoucan_ac = ""
                            glycosylation_type = obj["glycosylation_type"][0]
                            glycosylation_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                            amino_acid = obj["amino_acid"][0].lower()
                            amino_acid = amino_acid[0].upper() + amino_acid[1:]
                            combo_id = "%s,%s,%s,%s" %(pos,glytoucan_ac,glycosylation_type,amino_acid)
                            uniprot_url = url_dict["uniprotkb"] % (uniprotkb_ac + "#ptm_processing")
                            ev_obj = {"database":"UniProtKB", "id":id1, "url":uniprot_url}
                            ev_detail = "%s,%s,%s" % (combo_id,"UniProtKB", id1)
                            #If this combo has not been populated 
                            if combo_id not in glycosylation_dict:
                                glycosylation_dict[combo_id] = {
                                    "glytoucan_ac":glytoucan_ac
                                    ,"type":glycosylation_type
                                    ,"position":pos
                                    ,"residue":amino_acid
                                    ,"evidence":[{"database":"UniProtKB", "id":id1, "url":uniprot_url}]
                                }
                            elif ev_detail not in seen_glyco_ev:
                                #In case this combo has been populated, add UniProtKB evidence
                                glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                                seen_glyco_ev[ev_detail] = True


                            #Add pmid to evidence bucket
                            if "PubMed" in obj["source"]:
                                for ev in obj["evidence"]:
                                    pmid = ev
                                    pmid_url = url_dict["pubmed"] % (pmid)
                                    ev_obj = {"database":"PubMed", "id":pmid, "url":pmid_url}
                                    ev_detail = "%s,%s,%s" % (combo_id,"PubMed", pmid)
                                    if ev_detail not in seen_glyco_ev:
                                        glycosylation_dict[combo_id]["evidence"].append(ev_obj)
                                        seen_glyco_ev[ev_detail] = True


                #Add glycosylation from pdb
                main_key = "glycosylation_sites_pdb"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        cond_list = []
                        cond_list.append(len(obj["glycosylation_site"]) > 0)
                        cond_list.append(len(obj["amino_acid"]) > 0)
                        cond_list.append(len(obj["glycosylation_type"]) > 0)
                        if False not in cond_list:
                            pos = int(obj["glycosylation_site_uniprotkb"][0])
                            pdb_pos = int(obj["glycosylation_site"][0])
                            glytoucan_ac = ""
                            glycosylation_type = obj["glycosylation_type"][0]
                            glycosylation_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                            amino_acid = obj["amino_acid"][0].lower()
                            amino_acid = amino_acid[0].upper() + amino_acid[1:]
                            combo_id = "%s,%s,%s,%s" %(pos,glytoucan_ac,glycosylation_type,amino_acid)
                            pdb_id = obj["pdb_id"][0]
                            pdb_url = url_dict["pdb4glycosylation"] % (pdb_id)
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
                main_key = "pathway_accessions"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        for j in xrange(0, len(obj["kegg_id"])):
                            p_id = obj["kegg_id"][j]
                            p_name = ""
                            url = url_dict["kegg"] % (p_id)
                            pathway.append({"id":p_id, "name":p_name, "resource":"KEGG", "url":url})
                        for j in xrange(0, len(obj["reactome_id"])):
                            p_id = obj["reactome_id"][j]
                            p_name = obj["reactome_label"][j]
                            url = url_dict["reactome"] % (p_id)
                            pathway.append({"id":p_id, "name":p_name, "resource":"Reactome", "url":url})

            
                #Extract enzyme roles (glycosyltransferases, glycohydrolases)
                main_key = "glycosyltransferase"
                if main_key in db_obj[id1]:
                    keywords.append("glycosyltransferase-activity")
                main_key = "glycohydrolase"
                if main_key in db_obj[id1]:
                    keywords.append("glycohydrolase-activity")


                #Extract function
		function_list = []
                main_key = "function_annotation"
                dblabel_dict = {"refseq":"RefSeq", "uniprotkb":"UniProtKB"}
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        ann = obj["annotation"][0]
                        database_name = obj["database"][0].lower()
                        database_id = obj["database_id"][0]
                        database_url = url_dict[database_name] % (database_id)
                        if database_name == "uniprotkb":
                            database_id = database_id.split("-")[0]
                            database_url = url_dict[database_name] % (database_id + "#function")
                        database_label = dblabel_dict[database_name]
                        #ev_list = []
		        #for pmid in obj["evidence"]:
                        #    pmid_url = url_dict["pubmed"] % (pmid)
			#    ev_list.append({"database":"PubMed", "id":pmid, "url":pmid_url})
                        ev_list = [{"database":database_label, "id":database_id, "url":database_url}]
                        function_obj = { "annotation":ann, "url":database_url, "evidence":ev_list}
                        function_list.append(function_obj)


                #Extract disease
                doid2icd10 = {}
                disease_dict = {}
                main_key = "disease"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        do_id = obj["do_id"][0]
                        doid_name = obj["doid_name"][0]
                        doid_name = doid_name[0].upper() + doid_name[1:]
                        #doid_def = obj["doid_definition"][0]
                        icd10_id = obj["icd_10_cm_id"][0] if obj["icd_10_cm_id"] != [] else ""
                        doid2icd10[do_id] = icd10_id     
                        database_url = url_dict["do"] % (do_id)
                        if do_id not in disease_dict:
                            disease_dict[do_id] = {
                                "name":doid_name, "doid":do_id, "icd10":icd10_id, "url":database_url,
                                "evidence":[]
                            }
                        if obj["mim_id"] != []:
                            omim_id = obj["mim_id"][0]
                            ev_obj = {"database":"OMIM", "id":omim_id, "url":url_dict["omim"] % (omim_id)}
                            disease_dict[do_id]["evidence"].append(ev_obj)
                        disease_list.append(disease_obj)
                
                
                disease_list = []
                for do_id in disease_dict:
                    disease_list.append(disease_dict[do_id])





                #Extract mutation
                mutation_list = []
                main_key = "mutation"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        start_pos = int(obj["aa_pos"][0])
                        end_pos = start_pos
                        ref_aa = obj["ref_aa"][0]
                        alt_aa = obj["alt_aa"][0]
                        freq = obj["mut_freq"][0]
                        do_id = obj["do_id"][0]
                        doid_name = obj["do_name"][0].split("/")[1].strip()
                        doid_name = doid_name[0].upper() + doid_name[1:]
                        ann = "High frequency mutation in %s (DOID:%s) observed in %s subjects" % (doid_name,do_id,freq)
                        icd10_id = ""
                        database_url = url_dict["do"] % (do_id)
                        disease_obj = {"name":doid_name, "icd10":icd10_id, "doid":do_id, "url":database_url}
                        biomuta_url = url_dict["biomuta"] % (uniprotkb_ac)
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
                seen_ortholog = {}
                main_key = "ortholog"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        ortholog_canon = obj["uniprotkb_acc_canonical_ortholog"][0]
                        if ortholog_canon not in db_obj:
                            continue #ortholog is not canonical
                        ev_db_name = obj["evidence_dataset_name"][0]
                        ev_db_id = obj["evidence_dataset_id"][0]
                        ev_db_url = ""
                        if len(obj["evidence_dataset_id"]) == 1:
                            ev_db_url = url_dict[ev_db_name] % (obj["evidence_dataset_id"][0])
                        elif len(obj["evidence_dataset_id"]) == 2: 
                            ev_db_url = url_dict[ev_db_name] % (obj["evidence_dataset_id"][0],obj["evidence_dataset_id"][1])
                        if ortholog_canon not in seen_ortholog:
                            orthologs.append(
                                    {
                                        "uniprot_canonical_ac":ortholog_canon
                                        ,"protein_name":""
                                        ,"organism":obj["ortholog_species"][0]
                                        ,"tax_id":int(obj["ortholog_tax_id"][0])
                                        ,"evidence":[
                                            {
                                                "database":ev_db_name
                                                ,"id":ev_db_id
                                                ,"url":ev_db_url
                                            }
                                        ]
                                    }
                            )
                            seen_ortholog[ortholog_canon] = True
                        else:
                            ev_obj = {
                                    "database":ev_db_name
                                    ,"id":ev_db_id
                                    ,"url":ev_db_url
                            }
                            orthologs[-1]["evidence"].append(ev_obj)



                #Extract mutation
                isoform2locusobj = {}
                main_key = "genelocus"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        isoform = obj["uniprot_isoform_ac"][0]
                        isoform2locusobj[isoform] = {
                            "chromosome":obj["chr_id"][0]
                            ,"start_pos":int(obj["start_pos"][0])
                            ,"end_pos":int(obj["end_pos"][0])
                            ,"evidence":[
                                {
                                    "database":"Ensembl Transcript"
                                    ,"id":obj["transcript_id"][0]
                                    ,"url":url_dict["ensembl"] % (obj["transcript_id"][0])
                                }
                                ,{
                                    "database":"Ensembl Peptide"
                                    ,"id":obj["peptide_id"][0]
                                    ,"url":url_dict["ensembl"] % (obj["peptide_id"][0])
                                }
                            ]
                        }
                        

		#Extract from "crossref"
                crossref = []
                gene = []
                main_key = "uniprot_cross_references"
                bgee_id = ""
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        if len(obj["bgee_id"]) > 0:
                            bgee_id = obj["bgee_id"][0]
                        if "hgnc_id" in obj and "hgnc_label" in obj:
                            if len(obj["hgnc_id"]) > 0 and len(obj["hgnc_label"]) > 0:
                                hgnc_url = url_dict["hgnc"] % (obj["hgnc_id"][0])
                                gene.append({"name":obj["hgnc_label"][0],"url":hgnc_url})
                        if "mgi_id" in obj and "mgi_label" in obj:
                            if len(obj["mgi_id"]) > 0 and len(obj["mgi_label"]) > 0:
                                mgi_url = url_dict["mgi"] % (obj["mgi_id"][0])
                                gene.append({"name":obj["mgi_label"][0],"url":mgi_url})
                db_name, db_id, db_url = "UniProtKB", uniprotkb_ac, url_dict["uniprotkb"] % (uniprotkb_ac)
                crossref.append({"database":db_name, "id":db_id, "url":db_url})
                if refseq_ac != "":
                    db_name, db_id, db_url = "RefSeq", refseq_ac, url_dict["refseq"] % (refseq_ac)
                    crossref.append({"database":db_name, "id":db_id, "url":db_url})


                #Add pdb references
                main_key = "ac2pdb"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        for db_id in obj["pdb_id"]:
                            db_name = "PDB"
                            db_url = url_dict["pdb"] % (db_id)
                            crossref.append({"database":db_name, "id":db_id, "url":db_url})


                #Extract expression_normal
                seen_doid = {}
                seen_uberonid = {}
                expression_tissue_list = []
                main_key = "expression_normal"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        uberon_id = obj["uberon_anatomy_id"][0].split(":")[1]
                        if uberon_id in seen_uberonid:
                            continue
                        seen_uberonid[uberon_id] = True
                        uberon_name = obj["uberon_name"][0]
                        uberon_name = uberon_name[0].upper() + uberon_name[1:]
                        call = "yes" if obj["expression_call"][0] == "present" else "no"
                        url = url_dict["uberon"] % (uberon_id)
                        tissue_obj = {"name":uberon_name,"uberon":uberon_id, "url":url}
                        url = url_dict["bgee"] % (bgee_id)
                        ev_list = [{"database":"Bgee", "id":bgee_id, "url":url}]
                        expression_tissue_list.append({"tissue":tissue_obj, "present":call, "evidence":ev_list})
                        if uberon_id in uberonid2doidlist:
                            for do_id in uberonid2doidlist[uberon_id]:
                                seen_doid[do_id] = True

                #Extract expression_disease
                expression_disease_list = []
                main_key = "expression_disease"
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        do_id = obj["doid"][0]
                        do_name = obj["doname"][0]
                        parent_doid = obj["parent_doid"][0] if len(obj["parent_doid"]) > 0 else ""

                        #disregard if both doid and parent_doi are not mapped to uberonid
                        if do_id not in seen_doid and parent_doid not in seen_doid:
                            continue
                        icd10 = doid2icd10[do_id] if do_id in doid2icd10 else ""
                        significance = obj["significance"][0]
                        direction = obj["direction"][0]
                        url = url_dict["do"] % (do_id)
                        disease_obj = {"name":do_name,"doid":do_id, "icd10":icd10, "url":url}
                        url = url_dict["bioxpress"] % (uniprotkb_ac)
                        ev_list = [{"database":"BioXpress", "id":uniprotkb_ac, "url":url}]
                        expression_disease_list.append({"disease":disease_obj, 
                                "significant":significance, "trend":direction,"evidence":ev_list})
                                                                            




                #Extract sequence for canonical
                sequence = {
                    "sequence":seq_hash[id1]
                    ,"length":len(seq_hash[id1])
                }

                #Extract isoforms
                isoforms = []
                for obj in db_obj[id1]["idmapping"]:
                    for isoform in obj["reviewed_isoforms"] + obj["unreviewed_isoforms"]:
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
                            ,"url":url_dict["uniprot_isoform"] % (uniprotkb_ac, isoform)
                            ,"sequence":{
                                "sequence":isoform_seq
                                ,"length":len(isoform_seq)
                            },
                            "locus":locus_obj
                        })

               
                #Extract publication
                publication = []
                main_key = "citations"
                seen_publication = {}
                if main_key in db_obj[id1]:
                    for obj in db_obj[id1][main_key]:
                        pmid, title, journal_name, url = "", "", "", ""
                        if obj["pmid"] != []:
                            for j in xrange(0, len(obj["pmid"])):
                                pmid = obj["pmid"][j]
                                title = obj["title"][j] if j < len(obj["title"]) else ""
                                journal_name = obj["journal_name"][j] if j < len(obj["journal_name"]) else ""
                                url = url_dict["pubmed"] % (pmid)
                                if pmid not in seen_publication and "" not in [pmid, title, journal_name]:
                                    pub_obj = {"pmid":pmid, "title":title, "journal":journal_name,"url":url}
                                    publication.append(pub_obj)
                                    seen_publication[pmid] = True






                out_obj = {
                    "uniprot_canonical_ac":id1
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
                if "idmapping" in  db_obj[uniprot_canonical_ac]:
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



