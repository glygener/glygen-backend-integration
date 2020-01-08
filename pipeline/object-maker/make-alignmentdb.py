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

sys.path.append('../../glytools/')
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

    config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    species_obj = config_obj["speciesinfo"]




    data_dir = "reviewed/"
     

    #load workbook
    work_book = {}
    seq_hash = {}
    isoform2canon = {}
    canon2recname = {}
    file_list_obj = json.loads(open("../../conf/alignment_datasets.json", "r").read())
    #file_list_obj = json.loads(open("../../conf/toy_protein_datasets.json", "r").read())


    isoform2taxid = {}
    for species in file_list_obj:
        for sheet_name in file_list_obj[species]:
            in_file = data_dir + "/%s_protein_%s.csv" % (species, sheet_name)
            if sheet_name not in work_book:
                work_book[sheet_name] = {}
            work_book[sheet_name]["fields"] = []
            print "loading ", species, sheet_name
            libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "uniprotkb_canonical_ac")
    
            if sheet_name == "masterlist":
                for canon in work_book[sheet_name]["data"]:
                    isoform2taxid[canon] = species_obj[species]["taxid"]
                    for row in work_book[sheet_name]["data"][canon]:
                        for isoform in [row[-2], row[-1]]:
                            isoform2canon[isoform] = canon
                            isoform2taxid[isoform] = species_obj[species]["taxid"]

    sheet_name = "recnames"
    for canon in work_book[sheet_name]["data"]:
        for row in work_book[sheet_name]["data"][canon]:
            canon2recname[canon] = row[0]

    tax_info = {}
    sheet_name = "protein_homolog_clusters"
    work_book[sheet_name] = {}
    print "loading ", sheet_name
    in_file = data_dir + "/%s.csv" % (sheet_name)
    libgly.load_sheet_as_dict(work_book[sheet_name], in_file, ",", "uniprotkb_canonical_ac")
    for canon in work_book[sheet_name]["data"]:
        tmp_fl = work_book[sheet_name]["fields"]
        for tmp_row in work_book[sheet_name]["data"][canon]:
            tax_id = tmp_row[tmp_fl.index("tax_id")]
            tax_name = tmp_row[tmp_fl.index("tax_name")]
            tax_info[canon] = {"taxid": tax_id, "taxname":tax_name}




    data_grid = {}
    seen_row = {}
    for sheet_name in work_book:
        print "transforming", sheet_name
        for canon in work_book[sheet_name]["data"]:
            if canon not in data_grid:
                data_grid[canon] = {}
            if sheet_name not in data_grid[canon]:
                data_grid[canon][sheet_name] = []
            for row in work_book[sheet_name]["data"][canon]:
                s = canon + sheet_name + json.dumps(row)
                if s not in seen_row:
                    data_grid[canon][sheet_name].append(row)
                seen_row[s] = True


    
    aln_file_list = glob.glob("alignments/*.aln")


    members_dict = {}
    isoform2uniprotkbid = {}

    cls_count = 0
    for canon in data_grid:
            
        #Extract uniprotkb_id
        sheet_name = "info_uniprotkb"
        tmp_fl = work_book[sheet_name]["fields"]
        if sheet_name in data_grid[canon]:
            tmp_row = data_grid[canon][sheet_name][0]
            isoform2uniprotkbid[canon] = tmp_row[tmp_fl.index("uniprotkb_id")]
        
        sheet_name = "masterlist"
        tmp_fl = work_book[sheet_name]["fields"]
        if sheet_name in data_grid[canon]:
            for tmp_row in data_grid[canon][sheet_name]:
                isoform_one = tmp_row[tmp_fl.index("reviewed_isoforms")]
                isoform_two = tmp_row[tmp_fl.index("unreviewed_isoforms")]
                isoform2uniprotkbid[isoform_one] = isoform2uniprotkbid[canon]
                isoform2uniprotkbid[isoform_two] = isoform2uniprotkbid[canon]
        

        #Extract isoform clustering
        tmp_fl = work_book["masterlist"]["fields"]
        if "masterlist" in  data_grid[canon]:
            cls_id = "isoformset.uniprotkb.%s" % (canon)
            aln_file = "alignments/isoformset.uniprotkb.%s.aln" % (canon)
            if aln_file in aln_file_list:
                for row in data_grid[canon]["masterlist"]:
                    if cls_id not in members_dict:
                        members_dict[cls_id] = []
                        cls_count += 1
                    if row[tmp_fl.index("reviewed_isoforms")] != "":
                        members_dict[cls_id].append(row[tmp_fl.index("reviewed_isoforms")])
                    if row[tmp_fl.index("unreviewed_isoforms")] != "":
                        members_dict[cls_id].append(row[tmp_fl.index("unreviewed_isoforms")])
        
        #Extract homologset clusters      
        tmp_fl = work_book["protein_homolog_clusters"]["fields"]
        if "protein_homolog_clusters" in data_grid[canon]:
            tmp_fl = work_book["protein_homolog_clusters"]["fields"]
            for tmp_row in data_grid[canon]["protein_homolog_clusters"]:
                homolog_cluster_id = tmp_row[tmp_fl.index("homolog_cluster_id")]
                database = tmp_row[tmp_fl.index("database")]
                cls_id = "homologset.%s.%s" % (database,homolog_cluster_id)
                aln_file = "alignments/homologset.%s.%s.aln" % (database,homolog_cluster_id)
                if aln_file in aln_file_list:
                    if cls_id not in members_dict:
                        members_dict[cls_id] = []
                        cls_count += 1
                    members_dict[cls_id].append(canon)

    #Populate canon2cls
    canon2cls = {}
    for cls_id in members_dict:
        for canon in members_dict[cls_id]:
            if canon not in canon2cls:
                canon2cls[canon] = []
            canon2cls[canon].append(cls_id)



    ts = datetime.datetime.now(pytz.timezone('US/Eastern')).strftime('%Y-%m-%d %H:%M:%S %Z%z')
    out_obj_list = []
    for cls_id in members_dict:
        out_obj= { "cls_id":cls_id, "date":ts, "sequences":[],
            "algorithm":{
                "name":"clustalw2",
                "url":"http://www.clustal.org/clustal2/",
                "parameter":"Default Clustalw2 parameters"
            }
        }
        aln_file = "alignments/%s.aln" % (cls_id)
        aln_hash = load_msa(aln_file)
        for seq_ac in aln_hash:
            if seq_ac != "consensus":
                name = canon2recname[seq_ac] if seq_ac in canon2recname else ""
                seq_id, tax_id, tax_name = "", 0, ""
                if seq_ac in isoform2taxid:
                    tax_id = isoform2taxid[seq_ac]
                    tax_name = species_obj[str(tax_id)]["taxname"]
                if seq_ac in isoform2uniprotkbid:
                    seq_id = isoform2uniprotkbid[seq_ac]
                o = {"uniprot_ac":seq_ac, "uniprot_id":seq_id, "id":seq_ac, 
                        "name":name, "aln":aln_hash[seq_ac],
                        "tax_id":tax_id, "tax_name":tax_name}
                out_obj["sequences"].append(o)
        n1, n2, n3, pid = 0, 0, 0, "0.0%"
        if "consensus" in aln_hash:
            if aln_hash["consensus"].strip() != "":
                out_obj["consensus"] = aln_hash["consensus"]
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
        if record_count%1000 == 0:
            print " ... filtered %s objects" % (record_count)
    print " ... final created in: %s objects" % (record_count)

    #out_file = path_obj["jsondbpath"] + "alignmentdb.json"
    #with open(out_file, "w") as FW:
    #    FW.write("%s\n" % (json.dumps(fout_obj, indent=4)))
    #    print " ... final filtered in: %s objects" % (record_count)

    #fout_obj = {}
    for canon in canon2cls:
        obj = {"uniprot_canonical_ac":canon, "clusterlist":canon2cls[canon]}
        obj = order_obj(obj)
        #fout_obj[canon] = order_obj(obj)
        out_file = path_obj["jsondbpath"] + "clusterdb/%s.json" % (canon)
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(obj, indent=4)))



if __name__ == '__main__':
        main()


