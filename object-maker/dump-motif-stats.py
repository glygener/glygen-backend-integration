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



def main():


    file_list = glob.glob("jsondb/motifdb/*.json")
    gtc2mac, mac2gtc = {}, {}
    for json_file in file_list:
        doc = json.loads(open(json_file,"r").read())
        motif_ac, glytoucan_ac = doc["motif_ac"], doc["glytoucan_ac"]
        mac2gtc[motif_ac] = glytoucan_ac
        for obj in doc["glycans"]:
            gtc = obj["glytoucan_ac"]
            if motif_ac not in mac2gtc:
                mac2gtc[motif_ac] = {}
            mac2gtc[motif_ac][gtc] = True

    exit()


    parent_dict, canon_count, site_count = {}, {}, {}

    file_list = glob.glob("jsondb/sitedb/*.json")
    for json_file in file_list:
        doc = json.loads(open(json_file,"r").read())
        site_id = doc["id"]
        for obj in doc["glycosylation"]:
            glytoucan_ac = obj["glytoucan_ac"] if "glytoucan_ac" in obj else ""
            if glytoucan_ac not in site_count:
                site_count[glytoucan_ac] = {}
            site_count[glytoucan_ac][site_id] = True

    file_list = glob.glob("jsondb/glycandb/*.json")    
    for json_file in file_list:
        doc = json.loads(open(json_file,"r").read())
        glytoucan_ac = doc["glytoucan_ac"]
        for o in doc["glycoprotein"]:
            uniprot_canonical_ac = o["uniprot_canonical_ac"]
            if glytoucan_ac not in canon_count:
                canon_count[glytoucan_ac] = {}
            canon_count[glytoucan_ac][uniprot_canonical_ac] = True
        for o in doc["motifs"]:
            motif_ac = o["id"]
            if motif_ac not in parent_dict:
                parent_dict[motif_ac] = {}
            parent_dict[motif_ac][glytoucan_ac] = True
            if glytoucan_ac not in gtc2mac:
                gtc2mac[glytoucan_ac] = {}
            gtc2mac[glytoucan_ac][motif_ac] = True

    
    for glytoucan_ac in site_count:
        if glytoucan_ac == "":
            continue
        for site_id in site_count[glytoucan_ac]:
            motif_list = list(gtc2mac[glytoucan_ac].keys()) if glytoucan_ac in gtc2mac else [""]
            canon = site_id.split(".")[0]
            for motif_ac in motif_list:
                print ("%s,%s,%s,%s" % (glytoucan_ac, motif_ac, canon, site_id))
    exit()

    for glytoucan_ac in canon_count:
        for uniprot_canonical_ac in canon_count[glytoucan_ac]:
            motif_list = list(gtc2mac[glytoucan_ac].keys()) if glytoucan_ac in gtc2mac else [""]
            for motif_ac in motif_list:
                print ("%s,%s,%s" % (glytoucan_ac, motif_ac, uniprot_canonical_ac))

    exit()


    for motif_ac in parent_dict:
        motif_gtc = mac2gtc[motif_ac]
        tmp_list = []
        for glytoucan_ac in parent_dict[motif_ac]:
            n = 0
            if glytoucan_ac in canon_count:
                n =  len(list(canon_count[glytoucan_ac].keys()))
                tmp_list += list(canon_count[glytoucan_ac])
                for canon in canon_count[glytoucan_ac]:
                    print ("%s(%s),%s,%s" % (motif_ac, motif_gtc, glytoucan_ac,canon))
            #print ("%s(%s),%s,%s" % (motif_ac, motif_gtc, glytoucan_ac,n))
        total = len(list(set(tmp_list)))
        #print ("%s(%s),%s,%s" % (motif_ac, motif_gtc,"TOTAL", total))


if __name__ == '__main__':
    main()

