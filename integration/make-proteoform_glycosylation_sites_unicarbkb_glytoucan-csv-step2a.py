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
    in_file = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation.csv" % (species)
    load_dataframe(data_frame, "queryresults", in_file, ",")


    seen = {}

    time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')
    out_file1 = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation_types_%s.csv" % (species,time_stamp)
    out_file2 = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation_types_updated.csv" % (species)
    FW1 = open(out_file1, "w")
    FW2 = open(out_file2, "w")

    row = ["uckb_id", "glycosylation_type"]
    FW1.write("\"%s\"\n" % ("\",\"".join(row)))
    FW2.write("\"%s\"\n" % ("\",\"".join(row)))
    for ac in data_frame["queryresults"]:
        for row in data_frame["queryresults"][ac]:
            ac = ac.split("/")[-1]
            #uckb_id = row["uckb_id"][0]
            uckb_id = row["Id"][0]
            if uckb_id[0:5] == "comp_":
                continue
            if uckb_id in seen:
                continue
            url = "http://unicarbkb.org/structure/%s" % (uckb_id)
            html_file = "tmp/%s.html" % (uckb_id)
            cmd = "rm tmp/%s.html; wget %s -O %s" % (html_file, url,  html_file)
            x = commands.getoutput(cmd)
            if os.path.exists(html_file):
                cmd = "grep \"%s\" %s " % ("glycan structure</p>", html_file)
                glycan_type = commands.getoutput(cmd).strip().split(" ")[0][3:].lower()
                row = [uckb_id, glycan_type]
                FW1.write("\"%s\"\n" % ("\",\"".join(row)))
                FW2.write("\"%s\"\n" % ("\",\"".join(row)))
                seen[uckb_id] = True

    FW1.close()
    FW2.close()






if __name__ == '__main__':
	main()
