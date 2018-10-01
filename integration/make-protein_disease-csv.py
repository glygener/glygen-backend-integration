import os,sys,glob
import json
import csv
from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--infile",action="store",dest="infile",help="disease n3 file")
    (options,args) = parser.parse_args()
    for file in ([options.infile]):
        if not (file):
            parser.print_help()
            sys.exit(0)
	
	in_file = options.infile
	config_json = json.loads(open("../conf/config-1.json", "r").read())
	out_file = config_json["pathinfo"]["intermediate"]+"/"+"human_protein_disease.csv"
	
	doid_database = {}
	doid_others = {}
	data = []
	do_id =""
	dod_id =""
	rowcount = 0
	with open(in_file, "r") as FR:
		for line in FR:
			row = line.split(' ')
			if "DOID_" in row[0]:
				do_id = row[0].split('/')[-1].split("_")[1].replace(">","")
			if do_id not in doid_database:
				doid_database[do_id] = {}
				doid_database[do_id]["mim_id"] = []
				doid_database[do_id]["doid_name"] = ""
				doid_database[do_id]["doid_defination"] = ""
				doid_database[do_id]["doid_alternative_name"] = []
				doid_database[do_id]["icd_10_cm_id"] = ""
				doid_database[do_id]["icd_9_cm_id"] = ""
				doid_database[do_id]["kegg_id"] = ""
				doid_database[do_id]["mesh_id"] = ""
				doid_database[do_id]["umls_id"] = []
			if "OMIM:" in line:
				row = line.split(' ')
				mim_id = row[2].split('"')[1].split(":")[1]
				doid_database[do_id]["mim_id"] += [mim_id]
			if "annotatedTarget" in line:
				row = line.split('"')
				doid_defination = row[1]
				doid_database[do_id]["doid_defination"] = doid_defination
			row = line.split(' ')
			if "rdf-schema#label" in row[1]:
				row = line.split('"')
				doid_name = row[1]
				doid_database[do_id]["doid_name"] = doid_name
			row = line.split(' ')
			if "hasExactSynonym" in row[1]:
				row = line.split('"')
				doid_altname = row[1]
				doid_database[do_id]["doid_alternative_name"] += [doid_altname]
			if 'ICD10CM:' in line:
				row = line.split('"')
				icd_10_cm_id = row[1].split(":")[1]
				doid_database[do_id]["icd_10_cm_id"] = icd_10_cm_id
			if 'ICD9CM:' in line:
				row = line.split('"')
				icd_9_cm_id = row[1].split(":")[1]
				doid_database[do_id]["icd_9_cm_id"] =  icd_9_cm_id
			if 'KEGG:' in line:
				row = line.split('"')
				kegg_id = row[1].split(":")[1]
				doid_database[do_id]["kegg_id"] = kegg_id
			if 'MESH:' in line:
				row = line.split('"')
				mesh_id = row[1].split(":")[1]
				doid_database[do_id]["mesh_id"] = mesh_id
			if 'UMLS_CUI:' in line:
				row = line.split('"')
                umls_id = row[1].split(":")[1]
                doid_database[do_id]["umls_id"] += [umls_id]

	rowcount = 0
 	header_list = ['do_id','mim_id','doid_name','doid_definiation','doid_alternative_name',"icd_10_cm_id","icd_9_cm_id","kegg_id","mesh_id","umls_id","url"]	
	with open(out_file, 'w') as csvfile:  
		writer = csv.writer(csvfile)
		writer.writerow(header_list)
		for do_id in doid_database:									
			row = ['"' + do_id + '"']
			if doid_database[do_id]["mim_id"] != "":
				row.append("|".join (doid_database[do_id]["mim_id"]))
			else:
				row.append("")
			if doid_database[do_id]["doid_name"] != "":
				row.append(doid_database[do_id]["doid_name"]) 	
			else:
				row.append("")
			if doid_database[do_id]["doid_defination"] != "":
				row.append(doid_database[do_id]["doid_defination"])
			else:
				row.append("")
			if doid_database[do_id]["doid_alternative_name"] != "":
				row.append("|".join (doid_database[do_id]["doid_alternative_name"]))
			else:
				row.append("")
			if doid_database[do_id]["icd_10_cm_id"] != "":
				row.append(doid_database[do_id]["icd_10_cm_id"])
			else:
				row.append("")
			if doid_database[do_id]["icd_9_cm_id"] != "":
				row.append(doid_database[do_id]["icd_9_cm_id"])
			else:
				row.append("")
			if doid_database[do_id]["kegg_id"] != "":
				row.append(doid_database[do_id]["kegg_id"])
			else:
				row.append("")
			if doid_database[do_id]["mesh_id"] != "":
				row.append(doid_database[do_id]["mesh_id"])
			else:
				row.append("")
			if doid_database[do_id]["umls_id"] != "":
				row.append("|".join(doid_database[do_id]["umls_id"]))
			else:
				row.append("")
			row += ["http://disease-ontology.org/term/DOID%3A"+do_id]
			writer.writerow(row)

if __name__ == '__main__':
        main()
