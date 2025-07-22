#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import subprocess
from optparse import OptionParser
import libgly
import csvutil
import datetime
import pytz
import hashlib

__version__="1.0"
__status__ = "Dev"





def load_bcoid2fname(fname2bcoid, bcoid2fname):

    file_list = glob.glob("jsondb/bcodb/*.json")
    for in_file in file_list:
        bco_doc = json.loads(open(in_file, "r").read())
        bco_id = bco_doc["object_id"].split("/")[-2]
        desc = bco_doc["provenance_domain"]["name"]
        if "io_domain" in bco_doc:
            for obj in bco_doc["io_domain"]["output_subdomain"]:
                if "uri" in obj:
                    if "filename" in obj["uri"]:
                        file_name = obj["uri"]["filename"]
                        if file_name.strip() != "" and file_name.find(".stat.csv") == -1:
                            fname2bcoid[file_name] = bco_id
                            if bco_id not in bcoid2fname:
                                bcoid2fname[bco_id] = {}
                            bcoid2fname[bco_id][file_name] = desc
    return 



def main():

    

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-v","--ver",action="store",dest="ver",help="2.0.2")

    (options,args) = parser.parse_args()
    for file in ([options.ver]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    
    ver = options.ver


    global wrk_dir
    global field_dict
    global io_dict 
    global generated_dir
    global log_file
 

    DEBUG = False
    #DEBUG = True

    generated_dir = "/data/projects/glygen/generated/"
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker"
    reviewed_dir = wrk_dir + "/reviewed"
    log_file = "logs/make-tablemakerdb.log"
    
    msg = "make-tablemakerdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")
    
    fname2bcoid, bcoid2fname = {}, {}
    load_bcoid2fname(fname2bcoid, bcoid2fname)

    api_base_url = "https://dsapi.tst.glygen.org/drs/objects/"
    file_base_url = "https://data.glygen.org/ln2data/releases/data/v-%s/reviewed/" % (ver)
    ts = datetime.datetime.now(pytz.timezone('US/Eastern'))
    ts = ts.strftime('%Y-%m-%dT%H:%M:%SZ')
  

    flag_desc_dict = json.load(open("generated/misc/filter_flags.json"))
    flag_desc_dict["species_not_in_glygen"] = "Reported species is not in GlyGen"


    flag_dict = {}
    file_list = glob.glob("logs/*.global.log")
    for in_file in file_list:
        with open(in_file, "r") as FR:
            lcount = 0
            for line in FR:
                lcount += 1
                if lcount == 1:
                    f_list = line.split("\",\"")
                else:
                    row = line.split("\",\"")
                    for flag in row[-1].strip().split(";"):
                        flag = flag.replace("\"", "")
                        flag_dict[flag] = flag_desc_dict[flag] if flag in flag_desc_dict else ""


    species_obj = {}
    in_file = "generated/misc/species_info.csv"
    libgly.load_species_info(species_obj, in_file)
    species_dict = {}
    for tax_id in species_obj:
        if species_obj[tax_id]["is_reference"] == "yes":
            species_dict[str(tax_id)] = True


    extra_excluded_records = {}
    file_list = glob.glob("downloads/tablemaker/current/*")
    for in_file in file_list:
        file_name = in_file.split("/")[-1]
        df = {}
        csvutil.load_sheet(df, in_file, [], ",")
        f_list = df["fields"]
        src_row_idx = 0
        for row in df["data"]:
            src_row_idx += 1
            tax_id = row[f_list.index("Species")] if "Species" in f_list else ""
            if tax_id != "" and tax_id not in species_dict:
                if file_name not in extra_excluded_records:
                    extra_excluded_records[file_name] = []
                o = {"row_index":int(src_row_idx), "exclusion_flags":["species_not_in_glygen"]}
                extra_excluded_records[file_name].append(o)


    file_name_dict = {}
    file_list = glob.glob(reviewed_dir + "/*tablemaker*.csv")
    for in_file in file_list:
        if in_file.find("_citations_") != -1 or in_file.find(".stat.csv") != -1:
            continue
        file_name = in_file.split("/")[-1]
        df = {}
        csvutil.load_sheet(df, in_file, [], ",")
        f_list = df["fields"]
        for row in df["data"]:
            src_file_name = row[f_list.index("src_file_name")]
            if src_file_name not in file_name_dict:
                file_name_dict[src_file_name] = {}
            file_name_dict[src_file_name][file_name] = True
  

    seen_tablemaker_id = {} 
    for src_file_name in file_name_dict:
        tablemaker_id = src_file_name
        file_url = "https://data.glygen.org/ln2downloads/tablemaker/current/%s" % (src_file_name) 
        doc = {"tablemaker_id":tablemaker_id, "file_url":file_url,
            "usage":[], "exclusion_flag_desc":flag_desc_dict}
        for file_name in file_name_dict[src_file_name]:
            bco_id = "GLY_XXXXXX"
            dataset_url = "https://data.glygen.org/%s" % (bco_id)
            file_url = "https://data.glygen.org/ln2data/releases/data/v-%s/reviewed/%s" % (ver, file_name)
            obj = {
                "dataset_url":dataset_url, "dataset_file":file_name, "bco_id":bco_id, "file_url":file_url,
                "excluded_records":[]
            }         
            if file_name in extra_excluded_records:
                obj["excluded_records"] = extra_excluded_records[file_name]
            doc["usage"].append(obj)
            ds_log_file = "logs/%s.global.log" % (file_name.replace(".csv", ""))
            df = {}
            csvutil.load_sheet(df, ds_log_file, [], ",")
            f_list = df["fields"]
            for row in df["data"]:
                src_row_idx = row[f_list.index("src_row_idx")]
                flag_list = row[f_list.index("flag_list")].split(";")
                o = {"row_index":int(src_row_idx), "exclusion_flags":flag_list}
                obj["excluded_records"].append(o) 
     
        out_file = wrk_dir + "/tablemaker/%s.json" % (tablemaker_id)
        with open(out_file,"w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        seen_tablemaker_id[tablemaker_id] = True

        msg = "make-tablemakerdb: created file: jsondb/tablemakerdb/%s.json" % (tablemaker_id)
        csvutil.write_log_msg(log_file, msg, "a")

    doc = []
    for t_id in seen_tablemaker_id:
        doc.append("https://data.glygen.org/ln2data/releases/data/v-2.8.1/tablemaker/%s.json" %(t_id)) 
    out_file = wrk_dir + "/tablemaker/url_list.json"
    with open(out_file,"w") as FW:
        FW.write("%s\n" % (json.dumps(doc, indent=4)))
        


    cmd = "chmod -R 775 " + wrk_dir + "/jsondb/tablemakerdb/"
    x, y = subprocess.getstatusoutput(cmd)



if __name__ == '__main__':
    main()



