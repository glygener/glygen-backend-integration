import os
import csv
import sys
import json
import commands
import glob
import gzip

import libgly

from Bio import SeqIO




def main():

    tax_id, species = "7227", "fruitfly"
    in_file = "unreviewed/%s_protein_xref_refseq.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    refseq2canon = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        refseq = row[f_list.index("xref_id")]
        refseq2canon[refseq] = canon


    record_list = []
    file_list = glob.glob("downloads/ncbi/refseq/current/invertebrate.*.faa")
    for in_file in file_list:
        for record in SeqIO.parse(in_file, "fasta"):
            refseq_ac = record.id.split(".")[0]
            if refseq_ac in refseq2canon:
                record_list.append(record)

    out_file = "downloads/ncbi/refseq/current/refseq_protein_all_%s.faa" % (tax_id)
    with open(out_file, "w") as FW:
        SeqIO.write(record_list, FW, "fasta")


    record_list = []
    file_list = glob.glob("downloads/ncbi/refseq/current/invertebrate.*.gpff")
    for in_file in file_list:
        for record in SeqIO.parse(in_file, "genbank"):
            refseq_ac = record.id.split(".")[0]
            if refseq_ac in refseq2canon:
                record_list.append(record)

    out_file = "downloads/ncbi/refseq/current/refseq_protein_all_%s.gpff" % (tax_id)
    with open(out_file, "w") as FW:
        SeqIO.write(record_list, FW, "gb") 


if __name__ == '__main__':
        main()
