#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
from optparse import OptionParser
import commands




#######################################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-v","--ver",action="store",dest="ver",help="Dataversion [e.g: 1.11.2]")
    (options,args) = parser.parse_args()
    for key in ([options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    rel = options.ver

    file_list = glob.glob("/data/shared/glygen/releases/data/v-%s/jsondb/bcodb/*.json" % (rel))
    filename2bcoid = {}
    doc_dict = {}
    for bco_file in file_list:
        bco_id = bco_file.split("/")[-1].split(".")[0]
        rel = bco_file.split("/")[-4].split("-")[-1]
        doc = json.loads(open(bco_file, "r").read())
        if "io_domain" not in doc:
            continue
        if doc["io_domain"]["output_subdomain"] == []:
            continue
        file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"]
        file_name = file_name.strip()
        filename2bcoid[file_name] = bco_id



    data_dir = "/data/projects/glygen/generated/"
    status_dict = {"one":{}, "two":{}, "three":{}, "four":{}}
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
            if in_file in filename2bcoid:
                status_dict["four"][in_file] = True
            if obj["integration_status"]["status"] == "integrate_all":
                status_dict["three"][in_file] = True
        else:
            for species in obj["categories"]["species"]:
                in_file = "%s_%s_%s.%s" % (species,mol,ds_name, ds_format)
                status_dict["one"][in_file] = True
                if in_file in filename2bcoid:
                    status_dict["four"][in_file] = True
                if species not in obj["integration_status"]["excludelist"]:
                    status_dict["three"][in_file] = True

    file_list = glob.glob(data_dir + "datasets/reviewed/*.*")
    for in_file in file_list:
        if in_file.find("stat.") != -1 or in_file.find("GLY_") != -1:
            continue
        file_name = in_file.split("/")[-1]
        status_dict["two"][file_name] = True




    out_file = "/data/projects/glygen/generated/datasets/logs/dataset_status_v-%s.csv" % (rel)
    FW = open(out_file, "w")
    row = ["flag", "bco_id", "file_name"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))

    file_list = status_dict["one"].keys()
    file_list += status_dict["two"].keys()
    file_list += status_dict["three"].keys()
    file_list += status_dict["four"].keys()
    file_list = list(set(file_list))
    for in_file in file_list:
        flag = "1" if in_file in status_dict["one"] else "0"
        flag += "1" if in_file in status_dict["two"] else "0"
        flag += "1" if in_file in status_dict["three"] else "0"
        flag += "1" if in_file in status_dict["four"] else "0"
        bco_id = filename2bcoid[in_file] if in_file in filename2bcoid else "GLY_000000"
        row = [flag, bco_id, in_file]
        FW.write("\"%s\"\n" % ("\",\"".join(row)))
    FW.close()
    
    cmd = "chmod 775 " + out_file
    x = commands.getoutput(cmd)
    print "created logfile: %s" % (out_file)



if __name__ == '__main__':
    main()


