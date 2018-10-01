import os,sys,glob
import json
import csv
from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON
from Bio import SeqIO
from Bio.Seq import Seq

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-o","--organism",action="store",dest="organism",help="human or mouse")
    parser.add_option("-i","--infile",action="store",dest="infile",help="NT input file")
    
    (options,args) = parser.parse_args()
    for file in ([options.organism,options.infile]):
        if not (file):
            parser.print_help()
            sys.exit()

    organism = options.organism
    in_file = options.infile
    config_json = json.loads(open("../conf/config-1.json", "r").read())

    if organism == "human":
        fasta_file = config_json["pathinfo"]["fasta"]+"/human_protein_all.fasta"
    else:
        fasta_file = config_json["pathinfo"]["fasta"]+"/mouse_protein_all.fasta"

    out_file = config_json["pathinfo"]["intermediate"]+"/"+organism+"_protein_information.csv"

    seq_hash = {}
    with open(fasta_file, 'r') as fastafile:
        for record in SeqIO.parse(fastafile, "fasta"):
            canonical = record.id.split("|")[1]
            seq_hash[canonical] = record.seq.upper()

    ac2canonical = {}
    canonical = set()
    canonical_length = {}
    id_map_file = config_json["pathinfo"]["reviewed"]+"/"+ organism +"_protein_idmapping.csv"
    with open(id_map_file, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
        rowCount = 0
        for row in csvreader:
            rowCount += 1
            if rowCount == 1:
                continue
            ac = row[0].split("-")[0]
            ac2canonical[ac] = row[0]
            canonical_length[row[0]] = len(seq_hash[row[0]])
            canonical.add(row[0])

    canonical_uniprotkb_id = {}
    canonical_mass= {}

    with open(in_file, "r") as FR:
        for line in FR:
            row = line.strip().split(' ')
            if "http://purl.uniprot.org/core/mnemonic" in row[1]:
                ac = row[0].split("/")[-1].replace(">","")
                if ac in ac2canonical:
                    canonical_temp = ac2canonical[ac]
                    if canonical_temp not in canonical_uniprotkb_id:
                        canonical_uniprotkb_id[canonical_temp] = []
                    canonical_uniprotkb_id[canonical_temp].append(row[2].replace('"',""))
            if "<http://purl.uniprot.org/core/mass>" in row[1]:
                isoform = row[0].split("/")[-1].replace(">","")
                if isoform in canonical:
                    if isoform not in canonical_mass:
                        canonical_mass[isoform]=[]
		    canonical_mass[isoform].append(row[2].split("^^")[0].replace('"',""))

    headerList = ["uniprotkb_canonical_ac","uniprotkb_id","protein_mass","protein_length"]
    with open(out_file,'w') as outfile1:
        writer = csv.writer(outfile1, delimiter=',')
        writer.writerow(headerList)
        for item in canonical:
            if item in canonical_mass or item in canonical_uniprotkb_id:
                row=[item]
                if item in canonical_uniprotkb_id:
                    row+=["|".join(canonical_uniprotkb_id[item])]
                else:
                    row+=[""]
                if item in canonical_mass:
                    row+=["|".join(canonical_mass[item])]
                else:
                    row+=[""]
                row+=[canonical_length[item]]
                writer.writerow(row)

if __name__ == '__main__':
    main()
