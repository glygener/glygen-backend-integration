import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient

import csvutil




__version__="1.0"
__status__ = "Dev"



###############################
def main():


    in_dir = "/data/projects/glygen/generated/"

    in_file = in_dir + "misc/field_names.json"
    tmp_dict = json.loads(open(in_file, "r").read())

    seen = {}
    for f in tmp_dict:
        seen[f] = True

    ds_obj_list = json.loads(open(in_dir + "misc/dataset-masterlist.json", "r").read())     
   
    error_list = []
    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if ds_format != "csv":
            continue
        if obj["categories"]["species"] == []:
            file_name = "%s_%s.%s" % (mol, ds_name, ds_format)
            out_file = in_dir + "datasets/unreviewed/%s_%s.%s" % (mol, ds_name, ds_format)
            if os.path.isfile(out_file) == False:
                error_list.append({"err":"file does not exist", "file":file_name})
                continue
            cmd = "head -1 %s" % (out_file)
            f_list = commands.getoutput(cmd).strip().split("\",\"")
            for f in f_list:
                f = f.replace("\"", "")
                combo_id = "%s-%s" % (file_name, f)
                if f not in seen and combo_id not in seen:
                    o = {"err":"field not registered","file":file_name, "field":f}
                    error_list.append(o)
                    seen[combo_id] = True
        else:
            for species in obj["categories"]["species"]:
                file_name = "*_%s_%s.%s" % (mol,ds_name, ds_format)
                out_file = in_dir + "datasets/unreviewed/%s_%s_%s.%s" % (species,mol,ds_name, ds_format)
                if os.path.isfile(out_file) == False:
                    f_name = "%s_%s_%s.%s" % (species,mol,ds_name, ds_format)
                    o = {"err":"file does not exist","file":f_name}
                    error_list.append(o)
                    continue
                cmd = "head -1 %s" % (out_file)
                f_list = commands.getoutput(cmd).strip().split("\",\"")
                for f in f_list:
                    f = f.replace("\"", "")
                    combo_id = "%s-%s" % (file_name, f)
                    if f not in seen and combo_id not in seen:
                        o = {"err":"field not registered","file":file_name, 
                                "field":f}
                        error_list.append(o)
                        seen[combo_id] = True



    out_file = "/data/projects/glygen/generated/datasets/logs/field_names_qc.csv" 
    FW = open(out_file, "w")
    row = ["file_name", "field", "error"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))
    for obj in error_list:
        row = []
        for f in ["file", "field", "err"]:
            v = obj[f] if f in obj else ""
            row.append(v)
        FW.write("\"%s\"\n" % ("\",\"".join(row)))
    FW.close()   
    
    cmd = "chmod 775 " + out_file              
    x = commands.getoutput(cmd)                    
    print "created logfile: %s" % (out_file)    




if __name__ == '__main__':
	main()

