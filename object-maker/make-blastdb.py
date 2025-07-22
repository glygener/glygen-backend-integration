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
import subprocess


sys.path.append('../../glytools/')
import libgly
import csvutil


def main():

    global config_obj
    global path_obj
    global species_obj


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"


    species_obj = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)
    species_list = []
    for k in species_obj:
        obj = species_obj[k]
        if obj["short_name"] not in species_list and obj["is_reference"] == "yes":
            species_list.append(obj["short_name"])

    out_dir = "jsondb/blastdb"
    #out_dir = "generated/datasets/blastdb"
    #cmd = "rm -f %s*" % (out_dir)
    #x = subprocess.getoutput(cmd)

    prg = "/software/ncbi-blast-2.6.0+/bin/makeblastdb"
    for seqset in ["canonical", "all"]:
        for species in species_list:
            src = "reviewed/%s_protein_%ssequences.fasta" % (species,seqset)
            dst = out_dir + "/%ssequences_%s.fasta" % (seqset,species)
            out = out_dir + "/%ssequences_%s" % (seqset, species)
            cmd = "cp " + src + " " + dst
            x = subprocess.getoutput(cmd)
            cmd = "%s -in %s -out %s -parse_seqids -dbtype prot" % (prg, dst, out)
            x = subprocess.getoutput(cmd)
        cmd = "cat %s/%ssequences_*.fasta > %s/%ssequences_all.fasta" % (out_dir,seqset,out_dir,seqset)
        x = subprocess.getoutput(cmd)
        dst = out_dir + "/%ssequences_all.fasta" % (seqset)
        out = out_dir + "/%ssequences_all" % (seqset)
        cmd = "%s -in %s -out %s -parse_seqids -dbtype prot" % (prg, dst, out)
        x = subprocess.getoutput(cmd)

    log_file = "logs/make-blastdb.log"
    msg = "make-blastdb: created blast files"
    csvutil.write_log_msg(log_file, msg, "w")




if __name__ == '__main__':
    main()

