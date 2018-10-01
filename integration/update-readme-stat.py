import os,sys
import json
import csv
import glob
import commands

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

    file_list = glob.glob("reviewed/%s_*_*.csv" % (species))
    file_list += glob.glob("reviewed/%s_*_*.fasta" % (species))
    file_list += glob.glob("reviewed/%s_*_*.gp" % (species))
    file_list += glob.glob("reviewed/%s_*_*.aln" % (species))
    file_list += glob.glob("reviewed/%s_*_*.png" % (species))
    file_list += glob.glob("reviewed/%s_*_*.nt" % (species))

    for in_file in file_list:
        file_name = ".".join(os.path.basename(in_file).split(".")[:-1])
        file_ext = os.path.basename(in_file).split(".")[-1]


        prefix = "_".join(file_name.split(".")[0].split("_")[2:])
        stat_buffer = ""
        if file_ext == "csv":
            field_list = load_dataframe(data_frame, prefix, in_file, ",")
            value_list = {}
            main_id = field_list[0]
            for canon in data_frame[prefix]:
                for row_obj in data_frame[prefix][canon]:
                    if main_id not in value_list:
                        value_list[main_id] = []
                    value_list[main_id].append(canon)
                    for field in field_list[1:]:
                        for v in row_obj[field]:
                            if field not in value_list:
                                value_list[field] = []
                            value_list[field].append(v)
            stat_buffer = "\tStatistics:\n"
            for field in field_list:
                stat_buffer += "\t%10s : %s\n" % (len(sorted(set(value_list[field]))), field)
       

        readme_file_one = "reviewed/%s" % (file_name) + ".readme.manual.txt"
        readme_file_two = "intermediate/%s" % (file_name) + ".readme.txt"
        readme_file_three = "reviewed/%s" % (file_name) + ".readme.txt"


        manual_buffer = ""
        if os.path.exists(readme_file_one) == True:
            manual_buffer = open(readme_file_one, "r").read()
        
        with open(readme_file_two, "w") as FW:
            FW.write("%s\n" % (manual_buffer))
            FW.write("%s\n" % (stat_buffer))

        cmd = "/usr/bin/md5sum %s" % (readme_file_two)
        md5sum_one = commands.getoutput(cmd).strip().split(" ")[0]
        cmd = "/usr/bin/md5sum %s" % (readme_file_three)
        md5sum_two = commands.getoutput(cmd).strip().split(" ")[0]
        if md5sum_one != md5sum_two:
            cmd = "cp %s %s" % (readme_file_two, readme_file_three)
            x = commands.getoutput(cmd)
            print "updated ", file_name
            




if __name__ == '__main__':
	main()
