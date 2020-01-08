import os
import csv
import sys
import json
import commands
import glob
import gzip
import requests

sys.path.append('../../glytools/')
import libgly




def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]
 
    global path_obj
    path_obj = config_obj["pathinfo"]

    url_tmplt = "https://panelapp.genomicsengland.co.uk/api/v1/genes/?entity_name=%s"
        
    data_frame = {}
    in_file = "unreviewed/human_protein_xref_hgnc.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_name = row[f_list.index("database_label")]
        out_file = path_obj["downloads"] + "genomics_england/panels/%s.json" % (gene_name)
        if os.path.isfile(out_file) == True:
            continue
        url = url_tmplt % (gene_name)
        res = requests.get(url, verify=False)
        if res.content.strip() != "":
            res_obj = json.loads(res.content)
            if res_obj["results"] == []:
                continue
            with open(out_file, 'w') as FW:
                FW.write("%s\n" % (res.content))
                print "downloaded json for %s " % (gene_name)
        


if __name__ == '__main__':
        main()
