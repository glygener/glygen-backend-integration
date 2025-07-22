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

import libgly
import csvutil

import subprocess





def batch_records(file_list, record_type, record_info, log_file):

    batch_size = 1000000
    for in_file in file_list:
        file_name = in_file.split("/")[-1]
        doc = json.loads(open(in_file, "r").read())
        main_id_field = record_info[record_type]["mainid"]
        main_id = doc[main_id_field]
        main_doc = {}
        batch_dict = {}
        total = 0
        for k_one in doc:
            if type(doc[k_one]) is list:
                size = len(json.dumps(doc[k_one]))
                if size > batch_size:
                    tmp_list = []
                    n = len(doc[k_one])
                    batch_id = 1
                    for idx in range(0, len(doc[k_one])):
                        obj = doc[k_one][idx]
                        total += 1
                        tmp_list.append(obj)
                        if len(json.dumps(tmp_list)) > batch_size:
                            if batch_id not in batch_dict:
                                batch_dict[batch_id] = {
                                    "batchid":batch_id,
                                    "recordid":main_id,
                                    "recordtype":record_type,
                                    "sections":{}
                                }
                            batch_dict[batch_id]["sections"][k_one] = tmp_list
                            #print (k_one, batch_id, total)
                            tmp_list = []
                            batch_id += 1
                    if tmp_list != []:
                        if batch_id not in batch_dict:
                            batch_dict[batch_id] = {
                                "batchid":batch_id,
                                "recordid":main_id,
                                "recordtype":record_type,
                                "sections":{}
                            }
                        batch_dict[batch_id]["sections"][k_one] = tmp_list
                    main_doc[k_one] = []
                else:
                    main_doc[k_one] = doc[k_one]
            else:
                main_doc[k_one] = doc[k_one]
        out_file = "jsondb/%sdb/%s.json" % (record_type,main_id.replace("/", "_"))
        with open(out_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(main_doc, indent=4)))
        l = len(json.dumps(main_doc))
        for batch_id in batch_dict:
            batch_doc = batch_dict[batch_id]
            out_file = "jsondb/batchdb/%s.%s" % (record_type, file_name)
            out_file = out_file.replace(".json", "." + str(batch_id) + ".json")
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(batch_doc, indent=4)))
        msg = "make-batchdb: ... processed %s" % (in_file)
        csvutil.write_log_msg(log_file, msg, "a")
    
    return


       

