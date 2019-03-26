import os,sys
import json
import csv
import requests
import commands
import time
import datetime


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
    in_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_idmapping.csv" % (species)
    load_dataframe(data_frame, "proteinidmapping", in_file, ",")

    in_file = config_obj["pathinfo"]["downloads"] + "/biomuta/biomuta.csv"
    load_dataframe(data_frame, "mutation", in_file, ",")

    seen_canon = {}
    prefix = "proteinidmapping"
    for main_id in data_frame[prefix]:
        seen_canon[main_id] = True


    row = ["uniprot_canonical_ac","aa_pos","ref_aa","alt_aa","mut_freq","data_src","do_id","do_name"]
    print "\"%s\"" % ("\",\"".join(row))
    prefix = "mutation"
    for isoform in data_frame[prefix]:
        if isoform not in seen_canon:
            continue
        for obj in data_frame[prefix][isoform]:
            if int(obj["mut_freq"][0]) < 10:
                continue
            if obj["do_id"][0] in ["3963", "0070003", "3119"]:
                continue
            row = [isoform, obj["aa_pos"][0], obj["ref_aa"][0],obj["alt_aa"][0],
                    obj["mut_freq"][0],
                    obj["data_src"][0],obj["do_id"][0], obj["do_name"][0]]
            print "\"%s\"" % ("\",\"".join(row))




if __name__ == '__main__':
	main()
