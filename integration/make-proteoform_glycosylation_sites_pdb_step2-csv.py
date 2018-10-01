import os,sys
import json
import glob
import csv
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--infile",action="store",dest="infile",help="CSV input file")
    (options,args) = parser.parse_args()
    for file in ([options.infile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    in_file = options.infile
    file_name = in_file.split("/")[-1].split(".")[0].split("_new")[0] 
    config_json = json.loads(open("../conf/config-1.json", "r").read())
    if "human" in file_name:
        fasta_file = config_json["pathinfo"]["fasta"]+"/human_protein_all.fasta"
    else:
        fasta_file = config_json["pathinfo"]["fasta"]+"/mouse_protein_all.fasta"
    out_file = config_json["pathinfo"]["intermediate"]+"/"+file_name+".csv"
    log_file = config_json["pathinfo"]["intermediate"]+"/"+file_name+".log"

    seq_hash = {}
    with open(fasta_file, 'r') as fastafile:
        for record in SeqIO.parse(fastafile, "fasta"):
            canonical = record.id.split("|")[1]
            seq_hash[canonical] = record.seq.upper()

    fw2 = open(log_file, "w")
    rowcount = 0
    aa_convertor = {}
    with open(out_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        with open(in_file, 'r') as infile:
            csvreader = csv.reader(infile, delimiter=',', quotechar='|')
            for row in csvreader:
                rowcount+=1
                if rowcount == 1:
                    headerList = row
                    writer.writerow(headerList+["amino_acid_uniprotkb"])
                    fw2.write(",".join(headerList)+",amino_acid_from_seq,eror"+"\n")
                else:
                    if row[6]=="":
                        fw2.write(",".join(row)+"\n")
                    else:
                        if row[0] in seq_hash:
                            if len(seq_hash[row[0]])>= int(row[6]):
                                aa_uniprotkb = seq_hash[row[0]][int(row[6])-1]
                                if row[2]=="Asn" and aa_uniprotkb=="N":
                                    writer.writerow(row+["Asn"])
                                elif row[2]=="Ser" and aa_uniprotkb=="S":
                                    writer.writerow(row+["Ser"])
                                elif row[2]=="Thr" and aa_uniprotkb=="T":
                                    writer.writerow(row+["Thr"])
                                else:
                                    fw2.write(",".join(row)+","+seq_hash[row[0]][int(row[6])-1]+",there is mismatch"+"\n")
                            else:
                                fw2.write(",".join(row)+","+",glycosylation site is greater than length of protein sequence."+"\n")
                        else:
                            fw2.write(",".join(row)+","+"canonical accession not present in  the fasta file."+"\n")
    fw2.close()
 
if __name__ == '__main__':
        main()
