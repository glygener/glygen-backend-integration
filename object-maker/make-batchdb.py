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
import batchutil






def main():


    DEBUG = False
    #DEBUG = True

    log_file = "logs/make-batchdb.log"
    msg = "make-batchdb: started logging "
    csvutil.write_log_msg(log_file, msg, "w")

    batch_config_obj = json.loads(open("../conf/batch.json", "r").read())
    record_type_list = list(batch_config_obj.keys())
    if DEBUG:
        record_type_list = ["protein"]

    for record_type in record_type_list:
        file_list = glob.glob("jsondb/jumbodb/%sdb/*.json" % (record_type))
        batchutil.batch_records(file_list, record_type, batch_config_obj, log_file)

   

if __name__ == '__main__':
    main()

