import os,sys
import string
import commands 
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient





__version__="1.0"
__status__ = "Dev"




def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-v","--ver",action="store",dest="ver",help="Dataversion [e.g: 1.11.2]")
    (options,args) = parser.parse_args()
    for key in ([options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)


    doc_dict = {}
    for bco_file in glob.glob("/data/shared/glygen/releases/data/v-*/jsondb/corrected_bcodb/*"):
        bco_id = bco_file.split("/")[-1].split(".")[0]
        rel = bco_file.split("/")[-4].split("-")[-1]
        if bco_id not in doc_dict:
            doc_dict[bco_id] = {}
        if rel not in doc_dict[bco_id]:
            doc_dict[bco_id][rel] = json.loads(open(bco_file, "r").read())



    rel = options.ver

    logs_dir = "/data/projects/glygen/generated/datasets/logs/"
    out_file = logs_dir + "bco_io_v-%s.csv" % (rel)

    FW = open(out_file, "w")
    row = ["bco_id", "output_filename","input_filename","output_uri", "input_uri"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))

    for bco_id in doc_dict:
        if rel not in doc_dict[bco_id]:
            continue
        doc = doc_dict[bco_id][rel]
        out_obj_list = []
        for obj in doc["io_domain"]["output_subdomain"]:
            if obj["uri"]["filename"].find(".stat.csv") == -1:
                out_obj_list.append(obj["uri"])

        in_obj_list = []
        for obj in doc["io_domain"]["input_subdomain"]:
            in_obj_list.append(obj["uri"])
        if in_obj_list == []:
            in_obj_list = [{"filename":"not_specified", "uri":"not_specified"}]
        for i in range(len(out_obj_list)):
            for j in range(len(in_obj_list)):
                row = [bco_id, out_obj_list[i]["filename"], in_obj_list[j]["filename"], out_obj_list[i]["uri"], in_obj_list[j]["uri"]]
                FW.write("\"%s\"\n" % ("\",\"".join(row)))

    FW.close()

    cmd = "chmod 775 " + out_file
    x = commands.getoutput(cmd)
    print "created logfile: %s" % (out_file)






if __name__ == '__main__':
	main()

