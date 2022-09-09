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
from Bio import SeqIO


__version__="1.0"
__status__ = "Dev"



def load_io_dict ():
    
    reviewed_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/reviewed/"
    compiled_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/compiled/"
    misc_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/misc/"
    root_dict = {
        "/data/shared/glygen/downloads/":"https://data.glygen.org/ln2downloads/",
        "/data/projects/glygen/generated/datasets/unreviewed/":reviewed_root,
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
                if path not in seen:
                    url = path
                    for s in root_dict:
                        url = url.replace(s, root_dict[s])
                    if file_name not in io_dict:
                        io_dict[file_name] = []
                    io_dict[file_name].append({"path":path, "url":url})
                    seen[path] = True

    return io_dict


def update_io_domain(doc):

    r_root = "https://data.glygen.org/ln2releases/data/v-x.x.x/reviewed/"
    r_path = "/data/projects/glygen/generated/datasets/reviewed/"

    bco_id = doc["object_id"].split("/")[-2]
    in_obj_list = []
    out_obj_list = []
    out_filename_list = []
    for obj in doc["io_domain"]["output_subdomain"]:
        out_file_name = obj["uri"]["uri"].split("/")[-1]
        obj["uri"]["filename"] = out_file_name
        ignore = False
        for k in [".stat.csv", ".log"]:
            if out_file_name.find(k) != -1:
                ignore = True
        if ignore == True:
            continue
        file_name = out_file_name

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




def validate_io(doc, reviewed_dir):
    
    dir_map = {
        "ln2downloads":"/data/shared/glygen/downloads/",
        "ln2releases":"/data/projects/glygen/generated/datasets/",
        "ln2wwwdata":"/data/projects/glygen/generated/datasets/"
    }

    bco_id = doc["object_id"].split("/")[-2].replace(".json","")

    error_list = []
    base_url = "https:/data.glygen.org"

    for obj in doc["io_domain"]["output_subdomain"]:
        obj["uri"]["uri"] = obj["uri"]["uri"].replace("//","/")
        file_name = obj["uri"]["uri"].split("/")[-1]
        if file_name.find(".stat.csv") != -1 or file_name[-4:] == ".log":
            continue
        path = "/data/projects/glygen/generated/datasets/reviewed/" + file_name
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
                path = "/data/projects/glygen/generated/misc/" + file_name
            if obj["uri"]["uri"].find("/source/") != -1:
                path = "/data/shared/glygen/downloads/old_source/" + file_name
                obj["uri"]["uri"] = base_url + "/ln2downloads/old_source/" + file_name
        elif in_type_one == "ln2downloads":
            path = path.replace(in_type_one, dir_map[in_type_one])
            path = path.replace("/v-x.x.x/", "/current/")
        path = path.replace("//","/")
        if os.path.isfile(path) == False:
            error_list.append("%s:%s input file not found" % (bco_id, path))
    return error_list


def transform_bco(bco_id, doc):

    # bco_id -> object_id
    ver = doc["provenance_domain"]["version"]
    doc["object_id"] = "https://biocomputeobject.org/%s/%s" % (bco_id, ver)
    if "bco_id" in doc:
        doc.pop("bco_id")


    # remove orcid        
    if "contributors" in doc["provenance_domain"]:
        for obj in doc["provenance_domain"]["contributors"]:
            if "orcid" in obj:
                if obj["orcid"].strip() == "":
                    obj.pop("orcid")

    # bco_spec_version -> spec_version
    doc["spec_version"] = ""
    if "bco_spec_version" in doc:
        doc["spec_version"] = doc["bco_spec_version"]
        doc.pop("bco_spec_version")

    #https://github.com/glygener/glygen-backend-integration

    # checksum  ->  etag  
    doc["etag"] = ""
    if "checksum" in doc:
        doc["etag"] = doc["checksum"]
        doc.pop("checksum")
        
    # remove error_domain if it is empty
    if "error_domain" in doc:
        e = {'empirical_error': {}, 'algorithmic_error': {}} 
        if doc["error_domain"] == e:
            doc.pop("error_domain")
                    
    # remove parametric_domain if it is empty
    if "parametric_domain" in doc:
        if doc["parametric_domain"] == []:
            doc.pop("parametric_domain")
                    
    # remove orcid sha1_chksum and access_time
    if "execution_domain" in doc:
        if "software_prerequisites" in doc["execution_domain"]:
            for o in doc["execution_domain"]["software_prerequisites"]:
                if "uri" in o:
                    for k in ["sha1_chksum", "access_time"]:
                        if k in o["uri"]:
                            o["uri"].pop(k)
        if "script" in doc["execution_domain"]:
            for o in doc["execution_domain"]["script"]:
                if "uri" in o:
                    for q in ["sha1_chksum", "access_time"]:
                        if q in o["uri"]:
                            o["uri"].pop(q)
    #correct io_domain          
    if "io_domain" in doc:
        for k in ["input_subdomain", "output_subdomain"]:
            for o in doc["io_domain"][k]:
                if "uri" in o:
                    for q in ["sha1_chksum", "access_time"]:
                        if q in o["uri"]:
                            o["uri"].pop(q)

        
    #correct description_domain
    if "description_domain" in doc:
        if "pipeline_steps" in doc["description_domain"]:
            for obj in doc["description_domain"]["pipeline_steps"]:
                for k in ["input_list", "output_list"]:
                    for o in obj[k]:
                        for q in ["sha1_chksum", "access_time"]:
                            if q in o:
                                o.pop(q)
                for o in obj["prerequisite"]:
                    if "uri" in o:
                        for q in ["sha1_chksum", "access_time"]:
                            if q in o["uri"]:
                                o["uri"].pop(q)
                   
    if "extension_domain" in doc:
        obj = {
            "extension_schema": "http://www.w3id.org/biocompute/extension_domain/1.1.0/dataset/dataset_extension.json",
            "dataset_extension": {
                "additional_license": {
                    "data_license": "https://creativecommons.org/licenses/by/4.0/",
                    "script_license": "https://www.gnu.org/licenses/gpl-3.0.en.html"
                }
            }
        }
        dc = doc["extension_domain"]["dataset_categories"] if "dataset_categories" in doc["extension_domain"] else []
        obj["dataset_extension"]["dataset_categories"] = dc
        doc["extension_domain"] = [obj]
        if doc["extension_domain"][0]["dataset_extension"]["dataset_categories"] == []:
            doc["extension_domain"][0]["dataset_extension"]["dataset_categories"] = [
                {"category_value": "Reviewed", "category_name": "status"}
            ]

        
    return






def load_bcos_from_editor(bco_dict, filename2bcoid):

    db_obj = json.loads(open("../conf/db.json", "r").read())

    try:
        client = pymongo.MongoClient('mongodb://localhost:27017',
            username=db_obj["mongodbuser"],
            password=db_obj["mongodbpassword"],
            authSource=db_obj["mongodbname"],
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        dbh = client[db_obj["mongodbname"]]
        coll = "c_bco"
        for doc in dbh[coll].find({}):
            if "bco_id" in doc:
                doc["bco_id"] = doc["bco_id"].replace("DSBCO_", "GLY_")
                bco_id = doc["bco_id"].split("/")[-1].replace(".json","")
                if "io_domain" not in doc:
                    continue
                if "output_subdomain" not in doc["io_domain"]:
                    continue
                retired_flag = "active"
                if "dataset_categories" in doc["extension_domain"]:
                    for o in doc["extension_domain"]["dataset_categories"]:
                        if "category_name" in o and "category_value" in o:
                            cat_name = o["category_name"].replace(" ", "_")
                            cat_value = o["category_value"].lower().strip()
                            if cat_name in ["status", "dataset_status"]:
                                if cat_value.find("retired") != -1:
                                    retired_flag = "retired"
               
                ####
                #for obj in doc["io_domain"]["output_subdomain"]:
                #    file_name = obj["uri"]["filename"].strip()
                #    print (bco_id, file_name, retired_flag)
                ######

                if retired_flag == "retired":
                    continue
                
                for obj in doc["io_domain"]["output_subdomain"]:
                    file_name = obj["uri"]["filename"].strip()
                    ignore = False
                    for k in [".stat.csv", "GLY_0"]:
                        if file_name.find(k) != -1:
                            ignore = True
                    if ignore == True or file_name == "":
                        continue
                    filename2bcoid[file_name] = bco_id
                    bco_dict[bco_id] = doc
    except pymongo.errors.ServerSelectionTimeoutError as err:
        print ({"error_list":[{"error_code": "open-connection-failed"}]})
    except pymongo.errors.OperationFailure as err:
        print ({"error_list":[{"error_code": "mongodb-auth-failed"}]})
    return 


def load_bcos_from_release(bco_dict, filename2bcoid, rel):

    rel_dir = "/data/shared/glygen/releases/data/"
    file_list = glob.glob(rel_dir + "v-%s/jsondb/bcodb/*" % (rel))
    for bco_file in file_list:
        bco_id = bco_file.split("/")[-1].split(".")[0]
        doc = json.loads(open(bco_file, "r").read())
        for obj in doc["io_domain"]["output_subdomain"]:
            file_name = obj["uri"]["filename"].strip()
            ignore = False
            for k in [".stat.csv", "GLY_0"]:
                if file_name.find(k) != -1:
                    ignore = True
            if ignore == True or file_name == "":
                continue
            filename2bcoid[file_name] = bco_id
        bco_dict[bco_id] = doc

    return


def load_bcos_from_misc_dir(bco_dict, filename2bcoid):

    file_list = glob.glob("/data/projects/glygen/generated/misc/bcodb/*.json")
    for bco_file in file_list:
        bco_id = bco_file.split("/")[-1].split(".")[0]
        doc = json.loads(open(bco_file, "r").read())
        for obj in doc["io_domain"]["output_subdomain"]:
            file_name = obj["uri"]["filename"].strip()
            ignore = False
            for k in [".stat.csv", "GLY_0"]:
                if file_name.find(k) != -1:
                    ignore = True
            if ignore == True or file_name == "" or file_name[-4:] == ".log":
                continue
            filename2bcoid[file_name] = bco_id
        bco_dict[bco_id] = doc

    return




def get_current_file_list():
    generated_dir = "/data/projects/glygen/generated/"
    file_list = glob.glob(generated_dir + "/datasets/reviewed/*")
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


def file_status(file_list, fname2bcoid_one, fname2bcoid_two):
    
    for out_file_name in file_list:
        if out_file_name in fname2bcoid_one:
            bco_id = fname2bcoid_one[out_file_name]
            print ("1", bco_id, out_file_name)
        elif out_file_name in fname2bcoid_two:
            bco_id = fname2bcoid_two[out_file_name]
            print ("2", bco_id, out_file_name)
        else:
            print ("x", "N/A", out_file_name)


    return

def main():

    global wrk_dir
    global field_dict
    global io_dict 
    global generated_dir
   

    generated_dir = "/data/projects/glygen/generated/"
    wrk_dir = "/home/rykahsay/glygen-backend-integration/object-maker"
    reviewed_dir = wrk_dir + "/reviewed/"

    in_file = wrk_dir + "/generated/misc/field_names.json"
    field_dict = json.loads(open(in_file,  "r").read())

    io_dict = load_io_dict()


    bco_list = []
    old_rel = "1.12.3"

    file_name_list_in_reviewed_dir = get_current_file_list()

    #BCOs in the old release
    bco_dict_one, fname2bcoid_one = {}, {}
    load_bcos_from_release(bco_dict_one, fname2bcoid_one, old_rel)
    
    #BCOs from from editor
    bco_dict_two, fname2bcoid_two = {}, {}
    load_bcos_from_editor(bco_dict_two, fname2bcoid_two)
    
    #BCOs from generated/misc/bcodb
    bco_dict_three, fname2bcoid_three = {}, {}
    load_bcos_from_misc_dir(bco_dict_three, fname2bcoid_three)

    #for f in list(fname2bcoid_three.keys()):
    #    print (f, f in file_name_list_in_reviewed_dir)
    #exit()
    
    #for f in file_name_list:
    #    flag = "%s%s%s"% (f in fname2bcoid_one, f in fname2bcoid_two,f in fname2bcoid_three)
    #    print (f, flag)
    #exit()



    file_name_list = file_name_list_in_reviewed_dir
    #file_name_list = list(fname2bcoid_three.keys())
    #file_name_list = ["fruitfly_protein_reactions_rhea.csv"]

    cmd = "rm -f " + wrk_dir + "/jsondb/bcodb/*"
    x, y = subprocess.getstatusoutput(cmd)

    bco_dict = {}
    failed_list = []
    n_mapped,n_unmapped = 0, 0
    for out_file_name in file_name_list:
        if out_file_name in fname2bcoid_one:
            bco_id = fname2bcoid_one[out_file_name]
            bco_dict[bco_id] = bco_dict_one[bco_id]
            n_mapped += 1
        elif out_file_name in fname2bcoid_three:
            bco_id = fname2bcoid_three[out_file_name]
            transform_bco(bco_id, bco_dict_three[bco_id])
            bco_dict[bco_id] = bco_dict_three[bco_id]
            n_mapped += 1
        elif out_file_name in fname2bcoid_two:
            bco_id = fname2bcoid_two[out_file_name]
            transform_bco(bco_id, bco_dict_two[bco_id])
            bco_dict[bco_id] = bco_dict_two[bco_id]
            n_mapped += 1
        else:
            failed_list.append(out_file_name)
            n_unmapped += 1

  

    io_error_list = []
    for bco_id in sorted(bco_dict):
        doc = bco_dict[bco_id]
        update_io_domain(doc)
        io_error_list += validate_io(doc, reviewed_dir)
        if "_id" in doc:
            doc.pop("_id")

        #Set version in URLs to v-x.x.x
        doc["object_id"] = "https://biocomputeobject.org/%s/v-x.x.x" % (bco_id)
        doc["provenance_domain"]["version"] = "v-x.x.x"
        for obj in doc["description_domain"]["pipeline_steps"]:
            for o in obj["input_list"] + obj["output_list"]:
                parts = o["uri"].split("/")
                if len(parts) > 3:
                    parts[-3] = "v-x.x.x"
                    o["uri"] = "/".join(parts)
        out_file = wrk_dir + "/jsondb/bcodb/%s.json" % (bco_id)
        with open(out_file,"w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
        print ("Created BCO file: jsondb/bcodb/%s.json" % (bco_id))

    cmd = "chmod -R 775 " + wrk_dir + "/jsondb/bcodb/"
    x, y = subprocess.getstatusoutput(cmd)

    print ("\n\tDataset files mapped to BCOs: %s" % (n_mapped))
    print ("\tDataset files unmapped to BCOs: %s" % (n_unmapped))
    print ("\tBCOs created: %s" % (len(list(bco_dict.keys()))))

    if failed_list != []:
        print ("\nFailed to create BCOs for the following dataset files")
        for out_file_name in failed_list:
            print ("\t",out_file_name)

    if io_error_list != []:
        print ("\nIO validation report")
        for err in io_error_list:
            print (err)

if __name__ == '__main__':
    main()



