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
import datetime
import pytz
import subprocess

import libgly





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



############################
def load_msa(aln_file):

    aln_hash = {}
    with open(aln_file, "r") as FR:
        lcount = 0
        prev_ac = ""
        aln_start = 0
        for line in FR:
            lcount += 1
            if lcount > 3:
                if lcount == 4:
                    first_aln = line.split(" ")[-1].strip()
                    aln_start = line.find(first_aln)
                ac = line[:aln_start].strip()
                ac = "consensus" if prev_ac not in ["", "consensus"] and ac == "" else ac
                aln = line[aln_start:-1]
                if ac not in [""]:
                    if ac not in aln_hash:
                        aln_hash[ac] = ""
                    aln_hash[ac] += aln
                prev_ac = ac

    return aln_hash



#######################################
def main():

    config_obj = json.loads(open("../conf/config.json", "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    species_obj = {}
    in_file = path_obj["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)

    species_list = []
    for k in species_obj:
        if k.isdigit() or species_obj[k]["is_reference"] == "no":
            continue
        species_list.append(k)




    data_dir = "reviewed/"
     

    #load workbook
    work_book = {}
    seq_hash = {}
    isoform2canon = {}
    canon2recname = {}
    canon2submittedname = {}
    isoform2taxid = {}
    isoform2taxname = {}
    is_canon = {}
    seqac2id = {}

    sheet_list = ["masterlist", "recnames", "submittednames","info_uniprotkb"]
    for sheet_name in sheet_list:
        for species in species_list:
            in_file = data_dir + "/%s_protein_%s.csv" % (species, sheet_name)
            if os.path.isfile(in_file) == False:
                continue
            print ("make-alignmenetdb:", in_file)

            sheet_obj = {}
            sheet_obj["fields"] = []
            libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "uniprotkb_canonical_ac")
            if sheet_name == "masterlist":
                for canon in sheet_obj["data"]:
                    is_canon[canon] = True
                    isoform2taxid[canon] = species_obj[species]["tax_id"]
                    isoform2taxname[canon] = species_obj[species]["long_name"]
                    for row in sheet_obj["data"][canon]:
                        for isoform in [row[-2], row[-1]]:
                            isoform2canon[isoform] = canon
                            isoform2taxid[isoform] = species_obj[species]["tax_id"]
                            isoform2taxname[isoform] = species_obj[species]["long_name"]
            elif sheet_name == "recnames":
                for canon in sheet_obj["data"]:
                    for row in sheet_obj["data"][canon]:
                        canon2recname[canon] = row[0]
            elif sheet_name == "submittednames":
                for canon in sheet_obj["data"]:
                    for row in sheet_obj["data"][canon]:
                        canon2submittedname[canon] = row[0]
            elif sheet_name == "info_uniprotkb":
                for canon in sheet_obj["data"]:
                    for row in sheet_obj["data"][canon]:
                        seqac2id[canon] = row[sheet_obj["fields"].index("uniprotkb_id")]


    sheet_name = "homolog_clusters"
    homolog_clusters_df = {}
    homolog_clusters_df["fields"] = []
    in_file = data_dir + "/protein_%s.csv" % (sheet_name)
    libgly.load_sheet_as_dict(homolog_clusters_df, in_file, ",", "uniprotkb_canonical_ac")            

    members_dict = {}
    
    aln_file_list = glob.glob("alignments/homologset/*.aln")
    aln_file_list += glob.glob("alignments/isoformset/*/*.aln")
    cls_count = 0
    clsid2species = {}
    cls_id_list = []
    for aln_file in aln_file_list:
        cls_id = aln_file.split("/")[-1].replace(".aln", "")
        if cls_id.find("isoformset") != -1:
            species = aln_file.split("/")[-2]
            clsid2species[cls_id] = species
        cls_id_list.append(cls_id)


    ts = datetime.datetime.now(pytz.timezone('US/Eastern')).strftime('%Y-%m-%d %H:%M:%S %Z%z')
    out_obj_list = []
    for cls_id in cls_id_list:
        out_obj= { "cls_id":cls_id, "date":ts, "sequences":[],
            "algorithm":{
                "name":"clustalw2",
                "url":"http://www.clustal.org/clustal2/",
                "parameter":"Default Clustalw2 parameters"
            }
        }
        aln_file = "alignments/homologset/%s.aln" % (cls_id)
        if cls_id.find("isoformset") != -1:
            aln_file = "alignments/isoformset/%s/%s.aln" % (clsid2species[cls_id],cls_id)

        aln_hash = load_msa(aln_file)
        for seq_ac in aln_hash:
            if seq_ac != "consensus":
                name = canon2recname[seq_ac] if seq_ac in canon2recname else ""
                seq_id, tax_id, tax_name = "", 0, ""
                if seq_ac in isoform2taxid:
                    tax_id = isoform2taxid[seq_ac]
                    tax_name = species_obj[str(tax_id)]["long_name"]
                if seq_ac in seqac2id:
                    seq_id = seqac2id[seq_ac]
                o = {"uniprot_ac":seq_ac, "uniprot_id":seq_id, "id":seq_ac, 
                        "name":name, "aln":aln_hash[seq_ac],
                        "tax_id":tax_id, "tax_name":tax_name}
                out_obj["sequences"].append(o)
        n1, n2, n3, pid = 0, 0, 0, "0.0%"
        if "consensus" in aln_hash:
            out_obj["consensus"] = aln_hash["consensus"]
            if aln_hash["consensus"].strip() != "":
                n1 = aln_hash["consensus"].count("*")
                n2 = aln_hash["consensus"].count(":")
                n3 = aln_hash["consensus"].count(".")
                n = len(aln_hash["consensus"])
                pid = "%s" % (round(100.0*float(n1)/float(n), 2)) + "%"
        out_obj["identical_positions"] = n1
        out_obj["similar_positions"] = n2 + n3
        out_obj["identity"] = pid
        out_obj_list.append(out_obj)


    #fout_obj = {}
    record_count = 0
    for obj in out_obj_list:
        cond_list= []
        if False in cond_list:
            continue
        #fout_obj[obj["cls_id"]] = order_obj(obj)
        out_file = path_obj["jsondbpath"] + "/alignmentdb/%s.json" % (obj["cls_id"])
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        record_count += 1 
    print ("make-alignmenetdb: ... final created in: %s objects" % (record_count))




if __name__ == '__main__':
        main()


