#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict


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

        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
	root_obj =  config_obj[config_obj["server"]]["rootinfo"]


	db_obj = {}


        #data_dir = path_obj["htmlpath"] + "/csv/"
        data_dir = "/data/projects/glygen/generated/datasets/reviewed/"

        file_list_obj = json.loads(open(data_dir + "/glycan_file_list.json", "r").read())
        file_list1, file_list2 = [], []
        for file_name_part in file_list_obj["glycan"]:
            file_list1.append(data_dir + "/human_glycan_" + file_name_part + ".csv")
            file_list1.append(data_dir + "/mouse_glycan_" + file_name_part+ ".csv")
        for file_name_part in file_list_obj["proteoform"]:
            file_list2.append(data_dir + "/human_proteoform_" + file_name_part + ".csv")
            file_list2.append(data_dir + "/mouse_proteoform_" + file_name_part + ".csv")



        species2taxid= {"human":9606, "mouse":10090}
        taxid2name = {}
        map_file = data_dir + "/taxid2name.csv"
        with open(map_file, 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                row_count += 1
                if row_count == 1:
                    continue
                taxid2name[row[0]] = row[1].strip()


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


        glycotag_dict = {}
        glycotag_file = data_dir + "/glycan_subtype_tags.csv"
        with open(glycotag_file, 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                if row ==[]:
                    continue
                row_count += 1
                if row_count == 1:
                    continue
                glycotag_dict[row[0]] = row[1].strip()

        canon2hgncid = {}
        with open(data_dir + "/canonical2hgncid.csv", 'r') as FR:
            csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in csv_grid:
                row_count += 1
                if row_count == 1:
                    continue
                canon2hgncid[row[0]] = row[1].strip()


        ac2recname = {}
        for in_file in glob.glob(data_dir + "*_protein_recnames.csv"):
            file_name = os.path.basename(in_file)
            with open(in_file, 'r') as FR:
                dataFrame = csv.reader(FR, delimiter=',', quotechar='"')
                rowCount = 0
                for row in dataFrame:
                    rowCount += 1
                    if rowCount == 1:
                        continue
                    ac2recname[row[0]] = row[1]

        FW = open("tmp/glycan-fields.txt", "w")
	glyid2speciesindex = {}
        for in_file in file_list1 + file_list2:
                file_name = os.path.basename(in_file)
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
                                        main_id_index = field_list.index("glytoucan_acc")
                                        id1 = row[main_id_index]
                                        #id1 = row[0]
                                        if id1 not in db_obj:
						db_obj[id1] = {}
					if prefix not in db_obj[id1]:
						db_obj[id1][prefix] = []
					row_obj = {}
					for j in xrange(0,len(field_list)):
						id2 = field_list[j]
                                                if id1 == id2:
                                                    continue
                                                if prefix in ["sequences"]:
                                                    row_obj[id2] = [] if row[j].strip() == "" else [row[j].replace("\"","")]
                                                else:
                                                    row_obj[id2] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                                        db_obj[id1][prefix].append(row_obj)
                print " ...... done! %s (tax_id=%s)" % (file_name, tax_id)

        FW.close()

	
        
        record_count = 0 
        out_obj_list = []
	for id1 in db_obj:
                if "properties" not in db_obj[id1]:
			continue
                if "glycan_mass" not in db_obj[id1]["properties"][0]:
			continue
                
		#Start obj record, and populate simple values
		out_obj = {"glytoucan_ac":id1}
		
                #Extract properties
                out_obj["mass"], out_obj["number_monosaccharides"] = -1, -1
                if len(db_obj[id1]["properties"][0]["glycan_mass"]) > 0:
                        out_obj["mass"] = round(float(db_obj[id1]["properties"][0]["glycan_mass"][0]), 2)
                if len(db_obj[id1]["properties"][0]["monosaccharides"]) > 0:
                    out_obj["number_monosaccharides"] = int(db_obj[id1]["properties"][0]["monosaccharides"][0])

                #Extract sequences
                out_obj["iupac_extended"], out_obj["wurcs"], out_obj["glycoct"] = "", "", ""
                if len(db_obj[id1]["sequences"][0]["iupac_extended"]) > 0:
			out_obj["iupac_extended"] = db_obj[id1]["sequences"][0]["iupac_extended"][0]
		if len(db_obj[id1]["sequences"][0]["wurcs"]) > 0:
			out_obj["wurcs"] = db_obj[id1]["sequences"][0]["wurcs"][0]
		if len(db_obj[id1]["sequences"][0]["glycoct"]) > 0:
			out_obj["glycoct"] = db_obj[id1]["sequences"][0]["glycoct"][0]
                                
                #Populate motif
                out_obj["motifs"] = []                
                seen_motif = {}
                for mobj in db_obj[id1]["motif"]:
                    for j in xrange(0, len(mobj["glytoucan_acc_motif"])):
                        motif_id = mobj["glytoucan_acc_motif"][j]
                        motif_name = mobj["motif_name"][j]
                        if motif_id not in seen_motif:
                            out_obj["motifs"].append({"id":motif_id, "name":motif_name})
                        seen_motif[motif_id] = True


		#Populate classification
		cobj_list = []
		seen = {}
                for cobj in db_obj[id1]["classification"]:
                    if "type" in cobj and "subtype" in cobj:
                        if len(cobj["type"]) != 0 and len(cobj["subtype"]) != 0:
                            s = json.dumps(cobj)
			    if s not in seen:
			    	cobj_list.append(cobj)
			    seen[s] = True

		type_combo_list = []
                for cobj in cobj_list:
                        for j in xrange(0, len(cobj["type"])):
                            glycan_type = cobj["type"][j].strip()
                            glycan_subtype = cobj["subtype"][j].strip()
                            type_combo_list.append("%s|%s" % (glycan_type, glycan_subtype))

                out_obj["classification"] = []
                type_combo_list = sorted(set(type_combo_list))
                if len(type_combo_list) == 1:
                    glycan_type, glycan_subtype = type_combo_list[0].split("|")
                    type_url = url_dict["glycan_type"] % (glycotag_dict[glycan_type])
                    o = {"type":{"name":glycan_type, "url":type_url}}
                    if glycan_subtype != "no subtype":
                        subtype_url = url_dict["glycan_type"] % (glycotag_dict[glycan_subtype])
                        o["subtype"] = {"name":glycan_subtype, "url":subtype_url}
                    out_obj["classification"].append(o)
                else:
                    for type_combo in type_combo_list:
                        glycan_type, glycan_subtype = type_combo.split("|")
                        if glycan_subtype == "no subtype":
                            continue
                        type_url = url_dict["glycan_type"] % (glycotag_dict[glycan_type])
                        subtype_url = url_dict["glycan_type"] % (glycotag_dict[glycan_subtype])
                        o = {
                            "type":{"name":glycan_type, "url":type_url}
                            ,"subtype":{"name":glycan_subtype, "url":subtype_url}
                        }
                        out_obj["classification"].append(o)


                species_dict = {}
		for obj in db_obj[id1]["idmapping"]:
                        for tax_id in obj["tax_id"]:
                            if tax_id not in ["9606", "10090"]:
                                continue
                            tax_name = taxid2name[tax_id] if tax_id in taxid2name else ""
                            species_dict[tax_id] = {"taxid":int(tax_id), "name":tax_name, "evidence":[]}
                            for db_id in obj["glycomedb_id"]:
				url = url_dict["glycomedb"] % (db_id)
                                o = {"database":"GlycomeDB", "id":db_id, "url":url}
                                species_dict[tax_id]["evidence"].append(o)
                            for db_id in obj["uckb_id"]:
                                url = url_dict["uckb_glycan"] % (db_id)
                                o = {"database":"UniCarbKB","id":db_id, "url":url}
                                species_dict[tax_id]["evidence"].append(o)
                            #for db_id in obj["glyco_id"]:
                            #    url = url_dict["glycodb"] % (db_id)
                            #    o = {"database":"GlycoDB", "id":db_id, "url":url}
                            #    the resource doesn't exist anymore
                            #    species_dict[tax_id]["evidence"].append(o)
                out_obj["species"] = []
                for tax_id in species_dict:
                    out_obj["species"].append(species_dict[tax_id])


                seen_glycan = {}
		out_obj["crossref"] = []
                for obj in db_obj[id1]["idmapping"]:
                    if id1 in seen_glycan:
                        continue
                    seen_glycan[id1] = True
                    for db_id in obj["glycomedb_id"]:
                        url = url_dict["glycomedb"] % (db_id)
                        o = {"database":"GlycomeDB", "id":db_id, "url":url}
			out_obj["crossref"].append(o)
                    for db_id in obj["uckb_id"]:
                        url = url_dict["uckb_glycan"] % (db_id)
                        o = {"database":"UniCarbKB","id":db_id, "url":url}
			out_obj["crossref"].append(o)

                    for db_id in obj["pbch_id"]:
                        url, db_name = "", ""
                        if db_id[0:3] == "CID":
                            url = url_dict["pubchem_compound"] % (db_id.replace("CID", ""))
                            db_name = "PubChem Compound"
                        elif db_id[0:3] == "SID":
                            url = url_dict["pubchem_substance"] % (db_id.replace("SID", ""))
                            db_name = "PubChem Substance"
                        o = {"database":db_name, "id":db_id, "url":url}
                        out_obj["crossref"].append(o)
			
                    #The resource doesn't exist anymore
                    #for db_id in obj["glyco_id"]:
                    #    url = url_dict["glycodb"] % (db_id)
                    #    o = {"database":"GlycoDB", "id":db_id, "url":url}
                    #    out_obj["crossref"].append(o)


                seen_evlist = {}
                glycosylation_dict = {}
		if "glycosylation_sites_unicarbkb_glytoucan" in db_obj[id1]:
			for obj in db_obj[id1]["glycosylation_sites_unicarbkb_glytoucan"]:
                                for pos in obj["glycosylation_site"]:
                                        combo_id = "%s,%s" % (obj["uniprotkb_acc_canonical"][0], pos)
                                        o = {"uniprot_canonical_ac":obj["uniprotkb_acc_canonical"][0], 
                                                "protein_name":ac2recname[obj["uniprotkb_acc_canonical"][0]],
                                                "position":pos}
                                        o["evidence"] = []
                                        if combo_id not in glycosylation_dict:
                                            glycosylation_dict[combo_id] = o
                                            seen_evlist[combo_id] = []
                                        for uckb_id in obj["uckb_id"]:
                                                url = url_dict["uckb_glycan"] % (uckb_id)
                                                if uckb_id not in seen_evlist[combo_id]:
                                		    glycosylation_dict[combo_id]["evidence"].append({
						        	"database":"UniCarbKB",
               						    "id":uckb_id,
							    "url":url
						    })
                                                seen_evlist[combo_id].append(uckb_id)
                
                out_obj["glycoprotein"] = []                            
                for combo_id in glycosylation_dict:
                    out_obj["glycoprotein"].append(glycosylation_dict[combo_id])


		out_obj["enzyme"] = []
                if "enzyme" in db_obj[id1]:
                        seen_ac = {}
                        for obj in db_obj[id1]["enzyme"]:
                                gene_name = obj["gene_symbol_enzyme"][0] if len(obj["gene_symbol_enzyme"]) > 0 else ""
                                ac = obj["uniprotkb_acc_canonical_enzyme"][0] if len(obj["uniprotkb_acc_canonical_enzyme"]) > 0 else ""
                   		#desc = obj["enzyme_name"][0] if len(obj["enzyme_name"]) > 0 else ""
                                desc = ac2recname[ac] if ac in ac2recname else ""
                                db_id = canon2hgncid[ac]
                                gene_url = ""
                                if obj["tax_id_enzyme"][0] == "9606":
                                    gene_url = url_dict["hgnc"] % (db_id)
                                elif obj["tax_id_enzyme"][0] == "10090":
                                    gene_url = url_dict["mgi"] % (db_id)
                                o = {
					"protein_name":desc, 
					"uniprot_canonical_ac":ac,
					"gene":gene_name,
					"gene_link":gene_url
				}
                                if ac not in seen_ac:
				    out_obj["enzyme"].append(o)
                                seen_ac[ac] = True

                out_obj_list.append(out_obj)
                record_count += 1
                if record_count%1000 == 0:
                    print " ... compiled %s objects" % (record_count)

        print " ... final compiled %s objects" % (record_count)

        record_count = 0
        fout_obj = {}
	for obj in out_obj_list:
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



