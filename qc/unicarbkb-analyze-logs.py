import os,sys
import json
import csv

from optparse import OptionParser

import commands
import glob
import csvutil


__version__="1.0"
__status__ = "Dev"



###############################
def main():

    logs_dir = "/data/projects/glygen/generated/datasets/logs/unicarbkb/"
    file_list = glob.glob("%s/*.log" % (logs_dir))

    count_dict = {}
    s_list = ["flag_list"]
    for in_file in file_list:
        file_name = in_file.split("/")[-1].replace(".log", "")
        species = in_file.split("/")[-1].split("_")[0]
        file_name = "_".join(in_file.split("/")[-1].split("_")[1:]).replace(".log", "")

        if file_name not in count_dict:
            count_dict[file_name] = {}
        if species not in count_dict[file_name]:
            count_dict[file_name][species] = {}
        seen_row = {}
        data_frame = {}
        csvutil.load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            row_str = json.dumps(row[1:])
            #if row_str in seen_row:
            #    continue
            #seen_row[row_str] = True
            tmp_list = []
            for f in s_list:
                if row[f_list.index(f)].strip() != "":
                    tmp_list.append(row[f_list.index(f)])

            for qc_call in ["qc_passed", "qc_failed"]:
                if qc_call not in count_dict[file_name][species]:
                    count_dict[file_name][species][qc_call] = {"total":0}
            qc_call = row[f_list.index("qc_call")]
            if tmp_list != []:
                combo = "|".join(tmp_list)
                if combo not in count_dict[file_name][species][qc_call]:
                    count_dict[file_name][species][qc_call][combo] = 0
                count_dict[file_name][species][qc_call][combo] += 1
            count_dict[file_name][species][qc_call]["total"] += 1

    #print json.dumps(count_dict, indent=4)
    #exit()
    out_file = logs_dir + "/summary_stat.csv"
    FW = open(out_file, "w")
    row = ["input_file_name", "species","qc_call", "qc_flags", "record_count"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))


    for file_name in count_dict:
        for species in count_dict[file_name]:
            for qc_status in count_dict[file_name][species]:
                for qc_flag in count_dict[file_name][species][qc_status]:
                    n = count_dict[file_name][species][qc_status][qc_flag]
                    row = [file_name, species,qc_status, qc_flag, str(n)]
                    FW.write("\"%s\"\n" % ("\",\"".join(row)))
    FW.close()
    cmd = "chmod 775 " + out_file
    x = commands.getoutput(cmd)

















if __name__ == '__main__':
        main()


