#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import requests
import subprocess
from Bio import SeqIO

import csvutil


def main():

    global data_path
    global data_root
    global data_rel


    DEBUG = False
    #DEBUG = True

    bcoid2filename, filename2bcoid  = {}, {}
    file_list = glob.glob("jsondb/bcodb/*.json")
    for in_file in file_list:
        bco_id = in_file.split("/")[-1].split(".")[0]
        bco_doc = json.loads(open(in_file, "r").read())
        for obj in bco_doc["io_domain"]["output_subdomain"]:
            if "uri" in obj:
                if "filename" in obj["uri"]:
                    file_name = obj["uri"]["filename"]
                    if file_name.strip() != "" and file_name.find(".stat.csv") == -1:
                        filename2bcoid[file_name] = bco_id
                        if bco_id not in bcoid2filename:
                            bcoid2filename[bco_id] = {}
                        bcoid2filename[bco_id][file_name] = True


        
    log_file = "logs/make-recordsdb.log"
    msg = "make-recordsdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")
   
    file_count = 0
    record_count = 0
    bco_id_list = sorted(list(bcoid2filename.keys()))
    if DEBUG:
        #bco_id_list = json.loads(open("tmp/list.json", "r").read())
        tmp_list = []
        for bco_id in bco_id_list:
            idx = int(bco_id.split("_")[-1])
            if idx in [922]:
                tmp_list.append(bco_id)
        bco_id_list = tmp_list





    for bco_id in bco_id_list:
        idx = int(bco_id.split("_")[1])
        file_idx = 0
        if bco_id not in bcoid2filename:
            continue
        for file_name in bcoid2filename[bco_id]:
            in_file = "reviewed/%s" % (file_name)
            ext = file_name.split(".")[-1].lower()
            #if ext not in ["gz", "zip"]:
            #    continue
            #print (bco_id)
            if os.path.isfile(in_file) == False:
                continue
            row_idx = 0
            doc_list = []
            if ext in ["gz", "mp4", "gif"]:
                file_idx += 1
                file_count += 1
                row_idx += 1
                record_id = "%s_%s_%s" % (bco_id, file_idx, row_idx)
                row = ["file format %s (*.%s) file cannot be displayed as text" % (ext, ext)]
                doc_list.append({"recordid":record_id,"row":json.dumps(row)})
            elif ext in ["nt"]:
                file_idx += 1
                file_count += 1
                with open(in_file, "r") as FR:
                    lcount = 0
                    for line in FR:
                        row_idx += 1
                        record_id = "%s_%s_%s" % (bco_id, file_idx, row_idx)
                        row = [line.strip()]
                        doc_list.append({"recordid":record_id,"row":json.dumps(row)})
                        if row_idx%1000 == 0:
                            msg = "make-recordsdb: ... extracted %s objects from %s " % (row_idx, file_name)
                            csvutil.write_log_msg(log_file, msg, "a")
                        if lcount > 10000:
                            break
                        lcount += 1
            elif ext in ["fasta"]:
                file_idx += 1
                file_count += 1
                for record in SeqIO.parse(in_file, "fasta"):
                    row_idx += 1
                    record_id = "%s_%s_%s" % (bco_id, file_idx, row_idx)
                    row_str = json.dumps([record.description,str(record.seq.upper())])
                    doc_list.append({"recordid":record_id,"row":row_str})
                    if row_idx%1000 == 0:
                        msg = "make-recordsdb: ... extracted %s objects from %s " % (row_idx, file_name)
                        csvutil.write_log_msg(log_file, msg, "a")
            elif ext in ["csv", "tsv"]:
                file_idx += 1
                file_count += 1
                sep = "," if ext == "csv" else "\t"
                data_frame = {}
                csvutil.load_sheet(data_frame, in_file, [], sep)
                f_list = data_frame["fields"]
                for row in data_frame["data"]:
                    row_idx += 1
                    record_id = "%s_%s_%s" % (bco_id, file_idx, row_idx)
                    doc = {"recordid":record_id, "row":[]}
                    for f in f_list:
                        doc["row"].append(row[f_list.index(f)])
                    doc["row"] = json.dumps(doc["row"])
                    doc_list.append(doc)
                    if row_idx%1000 == 0:
                        msg = "make-recordsdb: ... extracted %s objects from %s " % (row_idx, file_name)
                        csvutil.write_log_msg(log_file, msg, "a")
            else:
                file_idx += 1
                file_count += 1
                row_idx += 1
                record_id = "%s_%s_%s" % (bco_id, file_idx, row_idx)
                row_str = json.dumps(["File format %s cannot be displayed, please download the file" %(ext)])
                doc_list.append({"recordid":record_id,"row":row_str})
                if row_idx%1000 == 0:
                    msg = "make-recordsdb: ... extracted %s objects from %s " % (row_idx, file_name)
                    csvutil.write_log_msg(log_file, msg, "a")

            if doc_list != []:
                out_dir = "jsondb/recordsdb/%s/" % (bco_id)
                if os.path.isdir(out_dir) == False:
                    cmd = "mkdir -p " + out_dir
                    x = subprocess.getoutput(cmd)
                grp = 0
                grp_dict = {}
                for idx in range(0, len(doc_list)):
                    doc  = doc_list[idx]
                    record_id = doc["recordid"]
                    if idx%1000 == 0:
                        grp += 1    
                        grp_dict[grp] = []
                    grp_dict[grp].append(doc)
                for grp in grp_dict:
                    out_file = out_dir + "%s_%s.json" % (bco_id, grp)
                    with open(out_file, "w") as FW:
                        FW.write("%s\n" % (json.dumps(grp_dict[grp], indent=4)))
                    record_count += 1
                msg = "make-recordsdb: ... created %s objects from %s " % (row_idx, file_name)
                csvutil.write_log_msg(log_file, msg, "a")
    
    msg = "make-recordsdb: ... finished total %s objects from %s files " % (record_count, file_count)
    csvutil.write_log_msg(log_file, msg, "a")


if __name__ == '__main__':
    main()



