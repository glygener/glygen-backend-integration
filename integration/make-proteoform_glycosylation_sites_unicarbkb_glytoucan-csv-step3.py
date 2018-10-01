import os,sys
import json
import csv
import requests
import commands
import time
import datetime

from Bio import SeqIO
from Bio.Seq import Seq

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


__version__="1.0"
__status__ = "Dev"


######################
def load_dataframe(data_frame, sheet_name, in_file, separator):

    data_frame[sheet_name] = {}
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='"')
        row_count = 0
        for row in csv_grid:
            row_count += 1
            if row_count == 1:
                field_list = row
            else:
                row_obj = {}
                for j in xrange(1,len(field_list)):
                    field_name  = field_list[j]
                    row_obj[field_name] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                main_id = row[0].strip()
                if main_id not in data_frame[sheet_name]:
                    data_frame[sheet_name][main_id] = []
                data_frame[sheet_name][main_id].append(row_obj)
    return



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


    data_frame = {}
    in_file = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation.csv" % (species)
    load_dataframe(data_frame, "glycosylation", in_file, ",")
    in_file = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation_types_updated.csv" % (species)
    load_dataframe(data_frame, "glycosylation_types", in_file, ",")
    
    in_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_idmapping.csv" % (species)
    load_dataframe(data_frame, "proteinidmapping", in_file, ",")
    in_file = config_obj["pathinfo"]["reviewed"] + "/%s_glycan_idmapping.csv" % (species)
    load_dataframe(data_frame, "glycanidmapping", in_file, ",")



    seq_hash = {}
    fasta_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_all.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())


    ac2canon = {}
    prefix = "proteinidmapping"
    for main_id in data_frame[prefix]:
        for obj in data_frame[prefix][main_id]:
            for isoform in obj["reviewed_isoforms"] + obj["unreviewed_isoforms"]:
                ac = isoform.split("-")[0]
                ac2canon[ac] = main_id

    uckbid2glytoucanac = {}
    prefix = "glycanidmapping"
    for main_id in data_frame[prefix]:
        for obj in data_frame[prefix][main_id]:
            for uckb_id in obj["uckb_id"]:
                uckbid2glytoucanac[uckb_id] = main_id

   
    uckbid2glycosylationtype = {}
    prefix = "glycosylation_types"
    for main_id in data_frame[prefix]:
        for obj in data_frame[prefix][main_id]:
            for t in obj["glycosylation_type"]:
                if t.find("n-linked") != -1 or t.strip() == "n":
                    if main_id not in uckbid2glycosylationtype:
                        uckbid2glycosylationtype[main_id] = {}
                    uckbid2glycosylationtype[main_id]["n-linked"] = True
                if t.find("o-linked") != -1:
                    if main_id not in uckbid2glycosylationtype:
                        uckbid2glycosylationtype[main_id] = {}
                    uckbid2glycosylationtype[main_id]["o-linked"] = True


    n_linked_aalist = ["Asn"]
    o_linked_aalist = ["Ser", "Thr"]

    aa_three2one = {"Asn":"N", "Ser":"S", "Thr":"T"}



    out_file = config_obj["pathinfo"]["intermediate"] + "/%s_proteoform_glycosylation_sites_unicarbkb_glytoucan.log" % (species)
    FW = open(out_file, "w")

    row = ["uniprotkb_acc_canonical","glycosylation_site","evidence",
            "uckb_id","glytoucan_acc","amino_acid","glycosylation_type"]
    print "\"%s\"" % ("\",\"".join(row))
    for ac in data_frame["glycosylation"]:
        for obj in data_frame["glycosylation"][ac]:
            ac = ac.split("/")[-1]
            #uckb_id = obj["uckb_id"][0]
            #amino_acid = obj["amino_acid"][0]
            #pos = obj["glycosylation_site"][0].split("^^")[0]
            #pm_id = obj["evidence"][0]
            
            #glytoucan_ac = uckbid2glytoucanac[uckb_id] if uckb_id in uckbid2glytoucanac else ""
                    
            uckb_id = obj["Id"][0]
            glytoucan_ac_new = obj["Toucan"][0].strip()
            glytoucan_ac = uckbid2glytoucanac[uckb_id].strip() if uckb_id in uckbid2glytoucanac else ""
            
            amino_acid = obj["TypeAminoAcid"][0]
            pm_id = obj["Pmid"][0]
            pos = obj["Position"][0].split("^^")[0]
            pos = int(pos) if pos != "" else 0
            canon = ""
            tv =  [False, False, False, False, False, False]
            if ac not in ac2canon:
                tv[0] = True
            else:
                canon = ac2canon[ac]
                if pos == 0 or pos >= len(seq_hash[canon]):
                    tv[1] = True
                elif aa_three2one[amino_acid] != seq_hash[canon][pos-1]:
                    tv[1] = True

            if uckb_id == "comp_":
                tv[2] = True
            if uckb_id not in uckbid2glycosylationtype:
                tv[3] = True
            
            tmp_list = []
            if uckb_id in uckbid2glycosylationtype:
                if "n-linked" in uckbid2glycosylationtype[uckb_id]:
                    tmp_list.append("n-linked")
                if "o-linked" in uckbid2glycosylationtype[uckb_id]:
                    tmp_list.append("o-linked")
                if tmp_list == ["n-linked"] and amino_acid not in n_linked_aalist:
                    tv[4] = True
                if tmp_list == ["o-linked"] and amino_acid not in o_linked_aalist:
                    tv[5] = True
            glycosylation_type = ";".join(tmp_list)

            row = [canon,str(pos),pm_id, uckb_id, glytoucan_ac_new,amino_acid,glycosylation_type]
            if True in tv:
                for j in xrange(0, len(tv)):
                    tv[j] = str(tv[j])
                row.append(";".join(tv))
                FW.write("\"%s\"\n" % ("\",\"".join(row)))
            else:
                print "\"%s\"" % ("\",\"".join(row))

    FW.close()

if __name__ == '__main__':
	main()
