#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import requests
import subprocess
import pymongo
from optparse import OptionParser
import libgly
import csvutil

from Bio import SeqIO


__version__="1.0"
__status__ = "Dev"



def load_io_dict ():
    
    reviewed_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/reviewed/"
    compiled_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/compiled/"
    misc_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/misc/"
    root_dict = {
        "/data/shared/glygen/downloads/":"https://data.glygen.org/ln2downloads/",
        "/data/shared/glyds/downloads/":"https://data.glygen.org/ln2downloads/",
        "/data/projects/glygen/generated/datasets/unreviewed/":reviewed_root,
        "/data/projects/glygen/generated/datasets/reviewed/":reviewed_root,
        "/data/projects/glygen/generated/datasets/compiled/":compiled_root,
        "/data/projects/glygen/generated/misc/":misc_root
    }








    io_dict = {}
    file_list = glob.glob(wrk_dir + "/usage/*.fu.log")
    for log_file in file_list:
        file_name = log_file.split("/")[-1].replace(".fu.log", "")
        with open(log_file, "r") as FR:
            seen = {}
            for line in FR:
                path = line.strip()
                path = path.replace("uniprot-proteome-sus_scrofa", "uniprot-proteome-sus-scrofa")
                path_parts = path.split("/")
                if path not in seen:
                    url = path
                    for s in root_dict:
                        url = url.replace(s, root_dict[s])
                    if path_parts[-2] in ["reviewed", "unreviewed"]:
                        url = reviewed_root + path_parts[-1]
                    if file_name not in io_dict:
                        io_dict[file_name] = []
                    io_dict[file_name].append({"path":path, "url":url})
                    seen[path] = True
    return io_dict


def update_dataset_categories(doc, bco_id):


    flag = False
    for f in bcoid2fname[bco_id]:
        if f.find("hcv1a_") != -1:
            flag = True
    if "extension_domain" not in doc:
        return
    
    new_obj_list = []
    seen = {}
    for obj in doc["extension_domain"]:
        if "dataset_extension" in obj:
            if "dataset_categories" in obj["dataset_extension"]:
                for o in obj["dataset_extension"]["dataset_categories"]:
                    if "category_name" in o and "category_value" in o:
                        if flag and o["category_name"] == "species" :
                            o["category_value"] = "Hepatitis C virus (isolate H77)"
                        s = "%s|%s" % (o["category_name"].lower(), o["category_value"].lower())
                        if s not in seen:
                            new_obj_list.append(o)
                            seen[s] = True
                obj["dataset_extension"]["dataset_categories"] = new_obj_list

    return 


def update_io_domain(doc):

    r_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/reviewed/"
    r_path = wrk_dir + "/reviewed/"

    ignore_file_list = [
        "chicken_proteoform_glycosylation_sites_unicarbkb.csv",
        "chicken_proteoform_glycosylation_sites_literature.csv",
        "chicken_proteoform_glycation_sites_uniprotkb.csv",
        "mouse_proteoform_citations_glycation_sites_uniprotkb.csv"
    ]


    bco_id = doc["object_id"].split("/")[-2]
    in_obj_list = []
    out_obj_list = []
    out_filename_list = []
    for obj in doc["io_domain"]["output_subdomain"]:
        #out_file_name = obj["uri"]["uri"].split("/")[-1]
        out_file_name = obj["uri"]["filename"].strip()
        if out_file_name == "":
            continue
        obj["uri"]["filename"] = out_file_name
        ignore = False
        for k in [".stat.csv", ".log"]:
            if out_file_name.find(k) != -1:
                ignore = True
        if ignore == True:
            continue
        file_name = out_file_name

        if "mediatype" not in obj:
            obj["mediatype"] = file_name.split(".")[-1]

        file_name = file_name.replace(".rdf.gz", ".gz")
        file_name = file_name.replace(".tar.gz", ".gz")
        file_name = ".".join(file_name.split(".")[:-1])
        
        if file_name not in out_filename_list:
            obj["uri"]["uri"] = r_root + out_file_name
            out_obj_list.append(obj)
            out_filename_list.append(file_name)
        input_count = 0
        if file_name in io_dict:
            for o in io_dict[file_name]:
                in_file_name = o["path"].split("/")[-1]
                in_url = o["url"]
                if in_url[0:6] == "/data/":
                    print ("\"%s\", bad in_url: %s" % (bco_id, in_url))
                    exit()
                if in_file_name not in ignore_file_list:
                    in_obj_list.append({"uri":{"filename":in_file_name, "uri":in_url}})
                    input_count += 1

    stat_out_obj_list = []
    for obj in out_obj_list:
        out_file_name = obj["uri"]["uri"].split("/")[-1]
        path = r_path + out_file_name
        file_ext = out_file_name.split(".")[-1]
        stat_rows = get_stat_rows(path)
        stat_file_name = ".".join(out_file_name.split(".")[:-1]) + ".stat.csv"
        stat_file_path = r_path + stat_file_name
        with open(stat_file_path, "w") as FW:
            row = ["unique_values","field_name", "field_description"]
            FW.write("\"%s\"\n" % ("\",\"".join(row)))
            for row in stat_rows:
                FW.write("\"%s\"\n" % ("\",\"".join(row)))
            stat_uri = r_root + stat_file_name
            o = {"mediatype": "csv", "uri":{"uri":stat_uri, "filename":stat_file_name}}
            stat_out_obj_list.append(o)


    #for o in in_obj_list:
    #    print (bco_id, o["uri"]["filename"])

    doc["io_domain"]["input_subdomain"] = in_obj_list
    doc["io_domain"]["output_subdomain"] = out_obj_list + stat_out_obj_list
    
    return




def get_stat_rows(in_file):

    file_ext = in_file.split(".")[-1]
    row_list = []
    if file_ext in ["csv", "tsv"]:
        sep = "," if file_ext == "csv" else "\t"
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, sep)
        stat_dict = {}
        for row in data_frame["data"]:
            #print (in_file)
            #print (row)
            #print (len(row), len(data_frame["fields"]))
            for j in range(0, len(row)):
                field = data_frame["fields"][j]
                if field not in stat_dict:
                    stat_dict[field] = {}
                stat_dict[field][row[j]] = True
        stat_obj = []
        for field_name in stat_dict:
            n = len(stat_dict[field_name].keys())
            d = field_dict[field_name]["description"] if field_name in field_dict else ""
            d = d.encode('ascii', 'ignore').decode('ascii')
            row = [str(n), field_name,d]
            row_list.append(row)
    elif file_ext == "fasta":
        n = len(list(SeqIO.parse(in_file, "fasta")))
        row = [str(n), "sequence","Amino acid sequence"]
        row_list.append(row)
    elif file_ext == "nt":
        n = 0
        with open(in_file, "r") as FR:
            for line in FR:
                n += 1
        row = [str(n), "triple","RDF triple"]
        row_list.append(row)
    else:
        row = ["1", "misc_field","Misc field type"]
        row_list.append(row)

    return row_list




def validate_io(doc):
    
    dir_map = {
        "ln2downloads":"/data/shared/glygen/downloads/",
        "ln2releases":wrk_dir + "/generated/datasets/",
        "ln2wwwdata": wrk_dir + "/generated/datasets/"
    }

    bco_id = doc["object_id"].split("/")[-2].replace(".json","")

    error_list = []
    base_url = "https:/data.glygen.org"

    for obj in doc["io_domain"]["output_subdomain"]:
        for p in  ["sha1_checksum", "access_time"]:
            if p in obj["uri"]:
                obj["uri"].pop(p)
        obj["uri"]["uri"] = obj["uri"]["uri"].replace("//","/")
        #file_name = obj["uri"]["uri"].split("/")[-1]
        file_name = obj["uri"]["filename"].strip()
        if file_name.find(".stat.csv") != -1 or file_name[-4:] == ".log":
            continue
        path = wrk_dir + "/reviewed/" + file_name
        if os.path.isfile(path) == False:
            error_list.append("%s:%s output file not found" % (bco_id, path))

    for obj in doc["io_domain"]["input_subdomain"]:
        obj["uri"]["uri"] = obj["uri"]["uri"].replace("//","/")
        obj["uri"]["uri"] = obj["uri"]["uri"].replace("ln2wwwdownloads","ln2downloads")
        obj["uri"]["uri"] = obj["uri"]["uri"].replace("gtc_gtc_02_04_2022","gtc_02_04_2022")
        obj["uri"]["uri"] = obj["uri"]["uri"].replace("doid2uberonid-mapping.csv", "doid2uberonid_mapping.csv")
        obj["uri"]["filename"] = obj["uri"]["filename"].replace("doid2uberonid-mapping.csv", "doid2uberonid_mapping.csv")
        file_name = obj["uri"]["uri"].split("/")[-1]
        in_type_one = obj["uri"]["uri"].split("/")[2]
        in_type_two = obj["uri"]["uri"].split("/")[-2]
        
        path = obj["uri"]["uri"]
        path = path.replace("https:/data.glygen.org","")
        path = path.replace("http:/data.glygen.org","")
        if in_type_one in ["ln2releases", "ln2wwwdata"]:
            path = dir_map[in_type_one]  + "/reviewed/" + file_name
            if obj["uri"]["uri"].find("/compiled/") != -1:
                path = dir_map[in_type_one]  + "/compiled/" + file_name
            if obj["uri"]["uri"].find("/misc/") != -1:
                path = wrk_dir + "/generated/misc/" + file_name
            if obj["uri"]["uri"].find("/source/") != -1:
                path = "/data/shared/glygen/downloads/old_source/" + file_name
                obj["uri"]["uri"] = base_url + "/ln2downloads/old_source/" + file_name
        elif in_type_one == "ln2downloads":
            path = path.replace(in_type_one, dir_map[in_type_one])
            path = path.replace("/v-x.x.x/", "/current/")
        path = path.replace("//","/")
        if os.path.isfile(path) == False and os.path.isdir(path) == False:
            error_list.append("%s:ERROR %s input file not found" % (bco_id, path))
    return error_list






def load_bcos_from_editor(bco_dict, filename2bcoid, bcoid2filename): 
    headers = {
        "Authorization":"Token 5467202e22e657dda9c19edb9a25936125163373",
        "Content-Type":"application/json"
    }
    req_obj = {
        "POST_api_objects_search":[
            {"type":"bco_id","search":"GLY_"}
        ]    
    }
    
    cmd =  'curl --request GET https://biocomputeobject.org/api/objects/search/'
    cmd += ' -H ' + '"Content-type:application/json"'
    cmd += ' -H ' + '"Authorization:Token 5467202e22e657dda9c19edb9a25936125163373"'
    cmd += ' -d ' + '\'{"POST_api_objects_search":[{"type":"bco_id","search":"GLY_"}]}\''
    cmd += ' -o logs/bco.txt'
    x = subprocess.getoutput(cmd)
    
    x = subprocess.getoutput("cat logs/bco.txt")
    bco_obj_list = json.loads(x)
    for bco_doc in bco_obj_list:
        bco_id = bco_doc["object_id"].split("/")[-2]
        if bco_id.find("GLY_") == -1:
            continue
        bco_ver = bco_doc["object_id"].split("/")[-1]
        if bco_ver == "DRAFT" and "contents" in bco_doc:
            status_list = []            
            if "extension_domain" in bco_doc["contents"]:
                for obj in bco_doc["contents"]["extension_domain"]:
                    if "dataset_extension" in obj:
                        o = obj["dataset_extension"]
                        if "dataset_categories" in o:
                            for oo in o["dataset_categories"]:
                                cat_name = oo["category_name"].lower() if "category_name" in oo else ""
                                cat_value = oo["category_value"].lower() if "category_value" in oo else ""
                                if cat_name == "status":
                                    status_list.append(cat_value)
            if "retired" in status_list:
                continue
            bco_dict[bco_id] = bco_doc["contents"]
            if "io_domain" in bco_doc["contents"]:
                for obj in bco_doc["contents"]["io_domain"]["output_subdomain"]:
                    if "uri" in obj:
                        if "filename" in obj["uri"]:
                            file_name = obj["uri"]["filename"]
                            if file_name.strip() != "" and file_name.find(".stat.csv") == -1:
                                filename2bcoid[file_name] = bco_id
                                if bco_id not in bcoid2filename:
                                    bcoid2filename[bco_id] = {}
                                bcoid2filename[bco_id][file_name] = True
        
    return 





def get_current_file_list():
    file_list = glob.glob(wrk_dir + "/reviewed/*")
    tmp_list = []
    for in_file in file_list:
        ignore = False
        for k in [".stat.csv", "GLY_0"]:
            if in_file.find(k) != -1:
                ignore = True
        if ignore == True:
            continue
        out_file_name = in_file.split("/")[-1]
        tmp_list.append(out_file_name)

    return tmp_list


def file_status(file_list, fname2bcoid_one, fname2bcoid):
    
    for out_file_name in file_list:
        if out_file_name in fname2bcoid_one:
            bco_id = fname2bcoid_one[out_file_name]
            print ("1", bco_id, out_file_name)
        elif out_file_name in fname2bcoid:
            bco_id = fname2bcoid[out_file_name]
            print ("2", bco_id, out_file_name)
        else:
            print ("x", "N/A", out_file_name)


    return

def get_bcoid2taxname(in_file):

    line_list = open(in_file, "r").read().split("\n")
    bcoid2taxname = {}
    for line in line_list[1:]:
        bco_id, tax_name = line.split(",")[0].replace("\"", ""), line.split(",")[-1].replace("\"", "")
        if bco_id not in bcoid2taxname:
            bcoid2taxname[bco_id] = []
        if tax_name not in bcoid2taxname[bco_id]:
            bcoid2taxname[bco_id].append({"category_name":"species", "category_value":tax_name})

    return bcoid2taxname


def get_bco_id_list(fname2bcoid, bcoid2fname):


    ds_file_status = {}
    file_name_list = get_current_file_list()
    for out_file_name in file_name_list:
        if out_file_name not in ds_file_status:
            ds_file_status[out_file_name] = []
        ds_file_status[out_file_name].append("in_fs")

    for bco_id in bcoid2fname:
        for out_file_name in bcoid2fname[bco_id]:
            if out_file_name not in ds_file_status:
                ds_file_status[out_file_name] = []
            ds_file_status[out_file_name].append("in_bco")
            


    bco_id_list = []
    for f in ds_file_status:
        bco_id = fname2bcoid[f] if f in fname2bcoid else "NO-BCO"
        flag_list = sorted(list(set(ds_file_status[f]))) if f in ds_file_status else []
        if flag_list != ['in_bco', 'in_fs']:
            msg = "make-bcodb: ERROR [%s, %s, %s]" % (f, bco_id, ";".join(flag_list))
            csvutil.write_log_msg(log_file, msg, "a")
        else:
            bco_id_list.append(bco_id)
    return bco_id_list



def main():

    global wrk_dir
    global field_dict
    global io_dict 
    global log_file
    global bcoid2fname


    #https://biocomputeobject.org/api/docs/#/BCO%20Management/api_objects_drafts_create_create

    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    in_file = wrk_dir + "/generated/misc/field_names.json"
    field_dict = json.loads(open(in_file,  "r").read())


    log_file = "logs/make-bcodb.log"
    msg = "make-bcodb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")
    

    io_dict = load_io_dict()
    #BCOs from from editor
    bco_dict, fname2bcoid, bcoid2fname = {}, {}, {}
    load_bcos_from_editor(bco_dict, fname2bcoid, bcoid2fname)

 
    DEBUG = False 
    #DEBUG = True

    in_file = wrk_dir + "/generated/misc/dataset2species.csv"
    bcoid2taxname = get_bcoid2taxname(in_file) 
    bco_id_list = get_bco_id_list(fname2bcoid, bcoid2fname) 
    if DEBUG:
        bco_id_list = json.loads(open("tmp/list.json", "r").read())
        bco_id_list = ["GLY_000716"]


    email_list = [
        "rsn13@gwu.edu",  "amandab2140@gwu.edu", "keeneyjg@gwu.edu",
        "mquartey@email.gwu.edu", "xavierh@email.gwu.edu"
    ]
    io_error_list = []
    for bco_id in sorted(bco_id_list):
        if bco_id not in bco_dict:
            continue
        idx = int(bco_id.split("_")[1])
        #if idx < 615:
        #    continue
        doc = bco_dict[bco_id]
        update_io_domain(doc)
        io_error_list += validate_io(doc)
        if "_id" in doc:
            doc.pop("_id")

        #Set version in URLs to v-x.x.x
        doc["object_id"] = "https://biocomputeobject.org/%s/v-x.x.x" % (bco_id)
        doc["provenance_domain"]["version"] = "v-x.x.x"
        if "description_domain" in doc:
            if "pipeline_steps" in doc["description_domain"]:
                for obj in doc["description_domain"]["pipeline_steps"]:
                    if "step_number" in obj:
                        obj["step_number"] = int(obj["step_number"])
                    obj["prerequisite"] = []
                    if "input_list" not in obj:
                        obj["input_list"] = []
                    if "output_list" not in obj:
                        obj["output_list"] = []
                    for o in obj["input_list"] + obj["output_list"]:
                        if o == None:
                            continue
                        if "uri" in o:
                            parts = o["uri"].split("/")
                            if len(parts) > 3:
                                parts[-3] = "v-x.x.x"
                                o["uri"] = "/".join(parts)
                    for k in ["input_list", "output_list"]:
                        for i in range(0, len(obj[k])):
                            if obj[k][i] == None:
                                continue
                            if obj[k][i] == {}:
                                obj[k][i]  = {"filename":"", "uri":""}
                            elif "uri" not in obj[k][i]:
                                obj[k][i]["uri"] = ""
                            for p in  ["sha1_checksum", "access_time"]:
                                if p in obj[k][i]:
                                    obj[k][i].pop(p)

        if "script" in doc["execution_domain"]:
            git_repo_url = "https://github.com/glygener/glygen-backend-integration/"
            for o in doc["execution_domain"]["script"]:
                fname = ""
                if "uri" in o["uri"]:
                    fname = o["uri"]["uri"].split("/")[-1]
                o["uri"]["uri"] = git_repo_url + "blob/master/dataset-maker/%s" % (fname)
        doc["execution_domain"]["software_prerequisites"] = [
            { "version": "2.7.5", "name": "Python",
                "uri": {"uri": "https://www.python.org/download/releases/2.7.5/"}
            }    
        ]
        doc["execution_domain"]["environment_variables"] = {}

        tax_name_map = {
            "sars coronavirus (sars-cov-1)":"Severe acute respiratory syndrome-related coronavirus",
            "sars coronavirus (sars-cov-2 or 2019-ncov)":"Severe acute respiratory syndrome coronavirus 2",
            "hepatitis c virus (genotype 1a, isolate h)":"Hepatitis c virus (isolate h77)",
            "hepatitis c virus (isolate h)":"Hepatitis c virus (isolate h77)",
            "hepatitis c virus (genotype 1b, isolate japanese)":"Hepatitis c virus (isolate japanese)",
            "saccharomyces cerevisiae (strain atcc 204508 / s288c)":"Saccharomyces cerevisiae S288C"
        }
        
        if "extension_domain" in doc:
            if len(doc["extension_domain"]) > 0:
                for obj in doc["extension_domain"]:
                    obj["extension_schema"] = obj["extension_schema"].replace("http://", "https://")
  
                obj = doc["extension_domain"][0]["dataset_extension"]
                if "dataset_categories" in obj:
                    o_list = []
                    for o in obj["dataset_categories"]:
                        if "category_name" in o and "category_value" in o:
                            if o["category_name"] == "species" and o["category_value"].strip() == "":
                                continue
                            if o["category_name"] == "protein":
                                o["category_name"] = "molecule"
                            if o["category_name"] == "tags":
                                o["category_name"] = "tag"
                            if o["category_name"] == "species" and o["category_value"].lower() in tax_name_map:
                                o["category_value"] = tax_name_map[o["category_value"].lower()]
                        o_list.append(o)
                    if bco_id in bcoid2taxname: 
                        o_list += bcoid2taxname[bco_id]
                    obj["dataset_categories"] = o_list

        cont_list = []
        for o in doc["provenance_domain"]["contributors"]:
            if "email" in o:
                if o["email"] in email_list:
                    continue
            cont_list.append(o)
        doc["provenance_domain"]["contributors"] = cont_list
       
        update_dataset_categories(doc, bco_id) 
        out_file = wrk_dir + "/jsondb/bcodb/%s.json" % (bco_id)
        with open(out_file,"w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        msg = "make-bcodb: created BCO file: jsondb/bcodb/%s.json" % (bco_id)
        csvutil.write_log_msg(log_file, msg, "a")

    cmd = "chmod -R 775 " + wrk_dir + "/jsondb/bcodb/"
    x, y = subprocess.getstatusoutput(cmd)


    if io_error_list != []:
        msg =  "make-bcodb: IO validation report"
        csvutil.write_log_msg(log_file, msg, "a")
        for err in io_error_list:
            msg = "make-bcodb: " + err
            csvutil.write_log_msg(log_file, msg, "a")

if __name__ == '__main__':
    main()



