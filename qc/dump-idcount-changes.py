import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
from Bio import SeqIO
import csvutil


__version__="1.0"
__status__ = "Dev"



def get_record_count(in_file):

    file_ext = in_file.split(".")[-1].lower()
    field_count, row_count, id_count = 1, 1, 1
    if file_ext in ["csv", "tsv"]:
        sep = "\t" if file_ext == "tsv" else ","
        field_count, row_count, id_count = csvutil.get_sheet_stats(in_file, sep)
    elif file_ext in ["fasta"]:
        id_count = len(list(SeqIO.parse(in_file, "fasta")))

    return field_count, row_count, id_count




###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-v","--ver",action="store",dest="ver",help="old reference ver")
    (options,args) = parser.parse_args()
    for key in ([options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)


    ver =  options.ver
    old_dir = "/data/shared/glygen/releases/data/v-%s/reviewed/" % (ver)
    new_dir = "/data/projects/glygen/generated/datasets/reviewed/"
    
    tmp_list = glob.glob(old_dir + "*.csv")
    tmp_list += glob.glob(old_dir + "*.tsv")
    tmp_list += glob.glob(old_dir + "*.fasta")
    old_file_list = []
    for f in tmp_list:
        if f.find("GLY_0") == -1 and f.find(".stat.csv") == -1:
            old_file_list.append(f)

    tmp_list = glob.glob(new_dir + "*.csv")
    tmp_list += glob.glob(new_dir + "*.tsv")
    tmp_list += glob.glob(new_dir + "*.fasta")
    new_file_list = []
    for f in tmp_list:
        if f.find("GLY_0") == -1 and f.find(".stat.csv") == -1:
            new_file_list.append(f)

    #old_file_list = sorted(old_file_list)[0:10]
    #new_file_list = sorted(new_file_list)[0:10]
    
    report_dict = {}
    for in_file in old_file_list:
        file_name = in_file.split("/")[-1]
        f_count, r_count, i_count = get_record_count(in_file)
        if file_name not in report_dict:
            report_dict[file_name] = {}
        report_dict[file_name]["old"] = i_count
    
    for in_file in new_file_list:
        file_name = in_file.split("/")[-1]
        f_count, r_count, i_count = get_record_count(in_file)
        if file_name not in report_dict:
            report_dict[file_name] = {}
        report_dict[file_name]["new"] = i_count
        

    out_file = "logs/idcount_changes.csv"
    FW = open(out_file, "w")
    row = ["idcount_change", "idcount_in_new", "idcount_in_old", "file_name"]
    FW.write("%s\n"  % (",".join(row)))
    for file_name in report_dict:
        n_old, n_new = 0, 0
        if "old" in report_dict[file_name]:
            n_old = report_dict[file_name]["old"]
        if "new" in report_dict[file_name]:
            n_new = report_dict[file_name]["new"]
        row = [str(abs(n_new - n_old)), str(n_new), str(n_old), file_name]
        FW.write("%s\n"  % (",".join(row)))
    FW.close()

    cmd = "chmod 775 %s" % (out_file)
    x = commands.getoutput(cmd)





if __name__ == '__main__':
	main()

