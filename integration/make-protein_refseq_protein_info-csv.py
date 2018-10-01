import os,sys
import util
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 
from Bio import SeqIO




__version__="1.0"
__status__ = "Dev"


###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--configfile",action="store",dest="configfile",help="config JSON file")
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

    (options,args) = parser.parse_args()
    for file in ([options.species, options.configfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    config_obj = json.loads(open(options.configfile, "r").read())
    species = options.species
    tax_id = "9606" if species == "human" else "10090"

    in_file = config_obj["pathinfo"]["downloads"] + "/refseq/protein/refseq_protein_all_%s.gp" % (tax_id)
    map_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_refseq_mapping.csv" % (species)

    refseq2canonical = {}
    with open(map_file, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_frame:
            row_count += 1
            refseq2canonical[row[3]] = row[0]

    data_frame = {}
    field = "xxx"
    with open(in_file, "r") as FR:
        for line in FR:
            newfield = line[0:12].strip() 
            field = newfield if len(newfield) > 0 else field
            value = line[12:].strip()
            if field == "VERSION":
                ac = line.split(" ")[-1].strip()
                data_frame[ac] = {"ac":ac, "summary":""}
                flag = False
            elif field == "COMMENT":
                if line.strip()[0:8] == "Summary:":
                    flag = True
                if line.strip() == "":
                    flag = False
                if flag:
                    data_frame[ac]["summary"] += line.strip() + " "


    
    row = ["uniprotkb_acc_canonical","p_refseq_acc_best_match","refseq_protein_name",
            "refseq_protein_length","refseq_protein_summary"]
    print "\"%s\"" % ("\",\"".join(row))
    for rec in SeqIO.parse(in_file, "genbank"):
        ac = rec.id
        for feat in rec.features:
            if feat.type in ["Protein"]:
                summary = data_frame[ac]["summary"]
                summary = summary.replace("\"", "`")
                summary = summary.replace("Summary: ", "")
                product = feat.qualifiers["product"][0]
                product = product.replace("\"", "`")
                seq_len = len(str(rec.seq))
                if ac in refseq2canonical:
                    uniprot_canonical_ac = refseq2canonical[ac]
                    row = [uniprot_canonical_ac, ac, product, str(seq_len), summary]
                    print "\"%s\"" % ("\",\"".join(row))



if __name__ == '__main__':
	main()


