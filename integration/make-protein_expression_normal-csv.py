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

    in_file = config_obj["pathinfo"]["downloads"] + "/bioxpress/%s_protein_expression_normal.csv" % (species)
    load_dataframe(data_frame, "normalexpression", in_file, ",")

    ac2canon = {}
    prefix = "proteinidmapping"
    for canon in data_frame[prefix]:
        ac = canon.split("-")[0]
        ac2canon[ac] = canon


    row = ["uniprot_canonical_ac","score","uberon_dev_id","sex",
            "uberon_anatomy_id","expression_call","uberon_name"]
    print "\"%s\"" % ("\",\"".join(row))
    prefix = "normalexpression"
    for ac in data_frame[prefix]:
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        for obj in data_frame[prefix][ac]:
            row = [canon, obj["expressionScore"][0], obj["uberonDevelopmentId"][0],obj["sex"][0],
                    obj["uberonAnatomyId"][0],
                    obj["expressionCall"][0],obj["uberon_name"][0]]
            print "\"%s\"" % ("\",\"".join(row))




if __name__ == '__main__':
	main()
