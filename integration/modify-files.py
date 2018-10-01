import os,sys
import json
import csv
import glob

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


######################
def load_dataframe(data_frame, sheet_name, in_file, separator):

    field_list = []
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
    return field_list





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
    tax_id = "9606" if species == "human" else "10090"



    data_frame = {}
    in_file = "reviewed/%s_blacklist.csv" % (species)
    field_list = load_dataframe(data_frame, "blacklist", in_file, ",")
        
    #file_list = glob.glob("reviewed/%s_*_*.csv" % (species))
    file_list = [
            "reviewed/mouse_protein_refseq_mapping.csv"
            ,"reviewed/human_protein_refseq_mapping.csv"
    ]


    for in_file in file_list:
        file_name = os.path.basename(in_file)
        prefix = "_".join(file_name.split(".")[0].split("_")[2:])
        field_list = load_dataframe(data_frame, prefix, in_file, ",")
        removed_count = 0
        out_file_one = "reviewed/%s" % (file_name)
        out_file_two = "logs/%s" % (file_name)
        FW1 = open(out_file_one, "w")
        FW2 = open(out_file_two, "w")
        FW1.write("\"%s\"\n" % ("\",\"".join(field_list)))
        for canon in data_frame[prefix]:
            for row_obj in data_frame[prefix][canon]:
                value_list = ["\"" + canon + "\""]
                for field in field_list[1:]:
                    tmp_list = []
                    for val in row_obj[field]:
                        tmp_list.append(val.strip())
                    value_list.append("\"" + "|".join(tmp_list) + "\"")
                row_line = ",".join(value_list)
                cond_list = [canon in data_frame["blacklist"]]
                for f in ["uniprotkb_acc", "uniprotkb_acc_canonical_enzyme"]:
                    if f in field_list:
                        for non_primary_canon in row_obj[f]:
                            cond_list.append(non_primary_canon in data_frame["blacklist"])
                if True in cond_list:
                    if removed_count == 0:
                        FW2.write("\"%s\"\n" % ("\",\"".join(field_list)))
                    removed_count += 1
                    FW2.write("%s\n" % (row_line))
                else:
                    FW1.write("%s\n" % (row_line))
        FW1.close()
        FW2.close()
        print removed_count, file_name





if __name__ == '__main__':
	main()
