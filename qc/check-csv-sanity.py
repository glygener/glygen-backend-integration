import os,sys
import json
import csv

import subprocess
import glob
import csvutil


###############################
def main():


    error_list = []
    in_dir = "/data/projects/glygen/generated/datasets/unreviewed/"
    file_list = glob.glob(in_dir+"/*.?sv")

    for in_file in file_list:
        flag = False
        for k in [".stat.", "GLY_"]:
            if in_file.find(k) != -1:
                flag = True
        if flag == True:
            continue
        if os.path.isfile(in_file) == True:
            delim = "," if in_file[-3:] == "csv" else "\t"
            data_frame = {}
            csvutil.load_sheet(data_frame, in_file, [], delim)
            f_list = data_frame["fields"]
            row_count = 0
            for row in data_frame["data"]:
                row_count += 1
                if len(f_list) != len(row):
                    err = "bad row: number of fields=%s, number of values=%s"
                    err = err % (len(f_list), len(row))
                    o = {"file":in_file, "row":row_count + 1, 
                            "error":err, "fields":f_list, "values":row}
                    error_list.append(o)

    if error_list == []:
        error_list = ["No errors found"]

    out_file = "/data/projects/glygen/generated/datasets/logs/sanity_qc.json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(error_list, indent=4)))
    
    cmd = "chmod 775 %s" % (out_file)
    x = subprocess.getoutput(cmd)

















if __name__ == '__main__':
        main()


