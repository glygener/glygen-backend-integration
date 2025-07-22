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
import csvutil




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
    in_file = "generated/misc/species_info.csv"
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

    log_file = "logs/make-clusterdb.log"
    msg = "make-clusterdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")

    sheet_list = ["masterlist", "recnames", "submittednames","info_uniprotkb"]
    for sheet_name in sheet_list:
        for species in species_list:
            in_file = data_dir + "/%s_protein_%s.csv" % (species, sheet_name)
            if os.path.isfile(in_file) == False:
                continue
            msg = "make-clusterdb: processing %s" % (in_file)
            csvutil.write_log_msg(log_file, msg, "a")

            sheet_obj = {}
            sheet_obj["fields"] = []
            libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "uniprotkb_canonical_ac")
            if sheet_name == "masterlist":
                for canon in sheet_obj["data"]:
                    is_canon[canon] = True
                    isoform2taxid[canon] = species_obj[species]["tax_id"]
                    isoform2taxname[canon] = species_obj[species]["common_name"]
                    for row in sheet_obj["data"][canon]:
                        for isoform in [row[-2], row[-1]]:
                            isoform2canon[isoform] = canon
                            isoform2taxid[isoform] = species_obj[species]["tax_id"]
                            isoform2taxname[isoform] = species_obj[species]["common_name"]
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
    for aln_file in aln_file_list:
        cls_id = aln_file.split("/")[-1].replace(".aln", "")
        if cls_id.find("isoformset") != -1:
            species = aln_file.split("/")[-2]
            clsid2species[cls_id] = species
            canon = cls_id.split(".")[-1]
            if canon not in is_canon:
                cmd = "rm -f alignments/isoformset/%s/isoformset.uniprotkb.%s.*" % (species,canon)
                x = subprocess.getoutput(cmd)
                msg = "REMOVED OLD FILES: alignments/isoformset/%s/isoformset.uniprotkb.%s.*" % (species,canon)
                csvutil.write_log_msg(log_file, msg, "a")
                continue
        aln_hash = load_msa(aln_file)
        if cls_id not in members_dict:
            members_dict[cls_id] = []
            if cls_count > 0 and cls_count%1000 == 0:
                msg = "make-clusterdb: processed %s clusters" % (cls_count)
                csvutil.write_log_msg(log_file, msg, "a")
            cls_count += 1
        for seq_ac in aln_hash.keys():
            if seq_ac != "consensus":
                members_dict[cls_id].append(seq_ac)
                if cls_id.find("isoform") != -1:
                    canon = isoform2canon[seq_ac]
                    seqac2id[seq_ac] = seqac2id[canon] 
    

    #Populate canon2cls
    canon2cls = {}
    for cls_id in members_dict:
        for seq_ac in members_dict[cls_id]:
            if seq_ac not in is_canon:
                continue
            if seq_ac not in canon2cls:
                canon2cls[seq_ac] = []
            canon2cls[seq_ac].append(cls_id)



    record_count = 0
    for canon in canon2cls:
        obj = {"uniprot_canonical_ac":canon, "clusterlist":canon2cls[canon]}
        obj = order_obj(obj)
        #fout_obj[canon] = order_obj(obj)
        out_file = "jsondb/clusterdb/%s.json" % (canon)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))
        record_count += 1


    msg = "make-clusterdb: ... final created in: %s objects" % (record_count)
    csvutil.write_log_msg(log_file, msg, "a")




if __name__ == '__main__':
        main()


