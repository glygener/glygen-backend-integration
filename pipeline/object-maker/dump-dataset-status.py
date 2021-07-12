#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob


sys.path.append('../../glytools/')
import libgly



#######################################
def main():


    data_dir = "/data/projects/glygen/generated/"
    status_dict = {"one":{}, "two":{}, "three":{}}
    ds_obj_list = json.loads(open(data_dir + "misc/dataset-masterlist.json", "r").read())
    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if ds_name in ["homolog_alignments", "isoform_alignments"]:
            continue
        if obj["categories"]["species"] == []:
            in_file = "%s_%s.%s" % (mol, ds_name, ds_format)
            status_dict["one"][in_file] = True  
            if obj["integration_status"]["status"] == "integrate_all":
                status_dict["three"][in_file] = True
        else:
            for species in obj["categories"]["species"]:
                in_file = "%s_%s_%s.%s" % (species,mol,ds_name, ds_format)
                status_dict["one"][in_file] = True
                if species not in obj["integration_status"]["excludelist"]:
                    status_dict["three"][in_file] = True

    file_list = glob.glob(data_dir + "datasets/reviewed/*.*")
    for in_file in file_list:
        if in_file.find("stat.") != -1:
            continue
        file_name = in_file.split("/")[-1]
        status_dict["two"][file_name] = True


    file_list = status_dict["one"].keys()
    file_list += status_dict["two"].keys()
    file_list += status_dict["three"].keys()
    file_list = list(set(file_list))
    for in_file in file_list:
        if in_file.find(".gz") != -1:
            continue
        flag = "1" if in_file in status_dict["one"] else "0"
        flag += "1" if in_file in status_dict["two"] else "0"
        flag += "1" if in_file in status_dict["three"] else "0"
        print flag, in_file


if __name__ == '__main__':
    main()


