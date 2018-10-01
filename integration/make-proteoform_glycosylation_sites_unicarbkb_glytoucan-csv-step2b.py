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
    
    in_file = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation_types_updated.csv" % (species)
    load_dataframe(data_frame, "glycosylation_types", in_file, ",")


    out_file = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation_types_updated.csv" % (species)
    FW = open(out_file, "w")
    row = ["uckb_id", "glycosylation_type"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))
    for uckb_id in data_frame["glycosylation_types"]:
        for obj in data_frame["glycosylation_types"][uckb_id]:
            if obj["glycosylation_type"] == []:
                url = "http://unicarbkb.org/structure/%s" % (uckb_id)
                html_file = "tmp/%s.html" % (uckb_id)
                cmd = "rm tmp/%s.html; wget %s -O %s" % (html_file, url,  html_file)
                x = commands.getoutput(cmd)
                if os.path.exists(html_file):
                    cmd = "grep \"%s\" %s " % ("glycan structure</p>", html_file)
                    glycan_type = commands.getoutput(cmd).strip().split(" ")[0][3:].lower()
                    if glycan_type != "":
                        obj["glycosylation_type"].append(glycan_type)
 
            glycan_type = obj["glycosylation_type"][0] if obj["glycosylation_type"] != [] else ""
            newrow = [uckb_id, glycan_type]
            FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
    FW.close()








if __name__ == '__main__':
	main()
