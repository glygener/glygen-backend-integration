import sys
import csv
import glob
import gzip


sys.path.append('../../glytools/')
import libgly



def main():
    
    do_dict = {}
    mondo_dict = {}
    in_file = "unreviewed/protein_disease_idmap.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        combo_id = "%s|%s" % (row[1], row[0])
        if combo_id not in do_dict:
            do_dict[combo_id] = []
        do_dict[combo_id].append(row[3])
        if combo_id not in mondo_dict:
            mondo_dict[combo_id] = []
        mondo_dict[combo_id].append(row[2])


    for combo_id in do_dict:
        print len(do_dict[combo_id]), combo_id



if __name__ == '__main__':
    main()




