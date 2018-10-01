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

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--infile",action="store",dest="infile",help="CSV input file")
    parser.add_option("-o","--organismType",action="store",dest="organism_type",help="The type of organism (human, mouse)")
    (options,args) = parser.parse_args()
    for file in ([options.infile,options.organism_type]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    organism = options.organism_type
    config_json = json.loads(open("../conf/config-1.json", "r").read())
    in_file1 = config_json["pathinfo"]["reviewed"]+"refseq_protein_sequence_10090.fasta"
    if "human" in organism:
        in_file1 = config_json["pathinfo"]["reviewed"]+"refseq_protein_sequence_9606.fasta"
        in_file2 = config_json["pathinfo"]["reviewed"]+"human_protein_canonical.fasta"
    else:
        in_file1 = config_json["pathinfo"]["reviewed"]+"refseq_protein_sequence_10090.fasta"
        in_file2 = config_json["pathinfo"]["reviewed"]+"mouse_protein_canonical.fasta"

    in_file3 = options.infile
    out_file = config_json["pathinfo"]["unreviewed"]+"/"+organism+"_protein_refseq_mapping.csv"

    seq_hash = {}
    for record in SeqIO.parse(in_file2, "fasta"):
        isoformAc = record.id.split("|")[1]
        seq_hash[isoformAc] = record.seq.upper()
    for record in SeqIO.parse(in_file1, "fasta"):
        refseqid = record.id.split(" ")[0].replace(">","")
        seq_hash[refseqid] = record.seq.upper()

    row_count = 0
    header_list = ["uniprotkb_acc_canonical","p_refseq_acc","n_refseq_acc","p_refseq_acc_best_match","n_refseq_acc_best_match"]
    with open(out_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header_list)
        with open(in_file3, 'r') as csvfile:
            csvReader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for row in csvReader:
                tmp_list = row[0:2]
                tmp_list += [row[2].replace('"',"")]
                row_count += 1
                if row_count == 1:
                    continue
                else:
                    seq_id1 = row[0].strip() 
                    if seq_id1 in seq_hash:
                        seq_str = ">" + seq_id1 + "\n" + seq_hash[seq_id1] + "\n"
                        in_file1 = "/tmp/acc_seq.fasta"
                        with open(in_file1, "w") as fw:
                            fw.write("%s" % (seq_str))
                        identity = []
                        refseq_accs = row[1].strip().split("|")
                        for seq_id2 in refseq_accs:
                            if seq_id2 in seq_hash:
                                seq_str = ">" + seq_id2 + "\n" + seq_hash[seq_id2] + "\n"
                                in_file2 = "/tmp/refseq_seq.fasta"
                                with open(in_file2, "w") as fw:
                                    fw.write("%s" % (seq_str))
                                outFile = "/tmp/outFile.txt"
                                cmd = "stretcher %s %s %s" % (in_file1, in_file2, outFile)
                                x = commands.getoutput(cmd)
                                with open(outFile, 'r') as txtFile:
                                    for line in txtFile:
                                        if "Identity" in line:
                                            identity += [float(line[line.find("(")+1:line.rfind(")")].replace("%",""))]
                        if identity==[]:
                            tmp_list += ["",""]
                        else:
                            maximum=0
                            for i,v in enumerate(identity):
                                if v>maximum:
                                    maximum=v
                                    index=i
                            tmp_list += [refseq_accs[index],row[2].split("|")[index].replace('"',"")]
                    else:
                        tmp_list += ["",""]
                    writer.writerow(tmp_list)

if __name__ == '__main__':
        main()
