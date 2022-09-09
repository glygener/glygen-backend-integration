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
import commands


sys.path.append('../../glytools/')
import libgly


def main():

    global config_obj
    global path_obj
    global species_obj


    config_file = "../conf/config.json"
    config_obj = json.loads(open(config_file, "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]

    data_dir = "reviewed/"


    species_obj = {}
    in_file = path_obj["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)
    species_list = []
    for k in species_obj:
        obj = species_obj[k]
        if obj["short_name"] not in species_list and obj["is_reference"] == "yes":
            species_list.append(obj["short_name"])

    cmd = "rm -f blastdb/*"
    x = commands.getoutput(cmd)

    prg = "/software/ncbi-blast-2.6.0+/bin/makeblastdb"
    for species in species_list:
        src = "reviewed/%s_protein_canonicalsequences.fasta" % (species)
        dst = "blastdb/canonicalsequences_%s.fasta" % (species)
        out = "blastdb/canonicalsequences_%s" % (species)
        cmd = "cp " + src + " " + dst
        x = commands.getoutput(cmd)
        
        cmd = "%s -in %s -out %s -parse_seqids -dbtype prot" % (prg, dst, out)
        x = commands.getoutput(cmd)

    cmd = "cat blastdb/canonicalsequences_*.fasta > blastdb/canonicalsequences_all.fasta"
    x = commands.getoutput(cmd)
    
    dst = "blastdb/canonicalsequences_all.fasta"
    out = "blastdb/canonicalsequences_all"
    cmd = "%s -in %s -out %s -parse_seqids -dbtype prot" % (prg, dst, out)
    x = commands.getoutput(cmd)




if __name__ == '__main__':
    main()

