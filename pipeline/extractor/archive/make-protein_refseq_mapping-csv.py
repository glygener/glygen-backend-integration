import os,sys
import string
import csv
import json
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment
import commands


sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"



def run_stretcher(seq_hash, seqid_one, seqid_two):

    cmd = "rm /tmp/seq_one.fasta /tmp/seq_two.fasta /tmp/stretcher.aln"
    x = commands.getoutput(cmd)

    with open("/tmp/seq_one.fasta", "w") as FW:
        FW.write("%s\n" % (">" + seqid_one + "\n" + seq_hash[seqid_one]))

    with open("/tmp/seq_two.fasta", "w") as FW:
        FW.write("%s\n" % (">" + seqid_two + "\n" + seq_hash[seqid_two]))

    cmd = "stretcher %s %s %s" % ("/tmp/seq_one.fasta", "/tmp/seq_two.fasta", "/tmp/stretcher.aln")
    x = commands.getoutput(cmd)
    
    cmd = 'grep "# Identity:" /tmp/stretcher.aln'
    x = commands.getoutput(cmd).strip().split("(")[-1].split("%")[0]
    return float(x)




    return





###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")


    (options,args) = parser.parse_args()
    for file in ([options.species]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    species = options.species
    species_info = {
        "human":{"taxid":9606, "taxname":"Homo sapiens"}
        ,"mouse":{"taxid":10090, "taxname":"Mus musculus"}    
    }

    config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())


    in_file1 = "reviewed/%s_protein_canonical.fasta" % (species)
    in_file2 = "reviewed/refseq_protein_sequence_%s.fasta" % (species_info[species]["taxid"])
    in_file3 = "downloads/ncbi/refseq/%s_protein_acc2refseqp.csv" % (species)
    out_file = "unreviewed/%s_protein_refseq_mapping.csv" % (species)


    ac2canon = {}
    seq_hash = {}
    for record in SeqIO.parse(in_file1, "fasta"):
        canon = record.id.split("|")[1]
        ac = canon.split("-")[0]
        seq_hash[canon] = record.seq.upper()
        ac2canon[ac] = canon

    for record in SeqIO.parse(in_file2, "fasta"):
        refseq_isoform = record.id.split(" ")[0].replace(">","")
        seq_hash[refseq_isoform] = record.seq.upper()

                            

    work_book = {}
    sheet_name = "refseqacmap"
    work_book[sheet_name] = {}
    libgly.load_sheet(work_book[sheet_name], in_file3, ",")


    row = ["uniprotkb_acc_canonical","refseq_acc","percent_identity"]
    print "\"%s\""  % ("\",\"".join(row))

    for row in work_book[sheet_name]["data"]:
        canon = ac2canon[row[0]]
        if canon in seq_hash:
            refseq_isoforms = row[1].split("|")
            identity_dict = {}
            for isoform in refseq_isoforms:
                if isoform in seq_hash:
                    identity_dict[isoform] = run_stretcher(seq_hash, canon, isoform)
            identity_dict = sorted(identity_dict.items(), key=lambda kv: kv[1], reverse=True)
            for k,v in identity_dict:
                row = [canon,k,v]
                print "\"%s\""  % ("\",\"".join(row))




if __name__ == '__main__':
        main()
