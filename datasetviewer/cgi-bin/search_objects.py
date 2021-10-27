#!/usr/bin/python
import os,sys
import string
import commands
import csv
import traceback
import cgi

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import util

import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"




###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    msg = "Input JSON text"
    parser.add_option("-j","--injson",action="store",dest="injson",help=msg)


    form_dict = cgi.FieldStorage()
    (options,args) = parser.parse_args()
    local_flag = False
    in_json = {}
    if len(form_dict.keys()) > 0:
        in_json = json.loads(form_dict["injson"].value) if "injson" in form_dict else {}
    else:
        local_flag = True
        for key in ([options.injson]):
            if not (key):
                parser.print_help()
                sys.exit(0)
        in_json = json.loads(options.injson)


    global config_json
    global db_obj
    global path_obj
    global client
    global data_path
    global data_root
    global current_release



    print "Content-Type: application/json"
    print   


    out_json = {}

    species_dict = {"human":"Homo sapiens", "mouse":"Mus musculus"}
    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        custom_config_json = json.loads(open("conf/config.custom.json", "r").read())
        db_obj = custom_config_json[config_json["server"]]["dbinfo"]
        path_obj = custom_config_json[config_json["server"]]["pathinfo"]
        root_obj = custom_config_json[config_json["server"]]["rootinfo"]
        current_release = "v-" + config_json["datarelease"]
        data_path = path_obj["htmlpath"] + "/ln2releases/"
        data_root = root_obj["htmlroot"] + "/ln2releases/"
    except Exception, e:
        print json.dumps({"taskstatus":0, "errormsg":"Loading config failed!"})
        sys.exit()

    try:
        client = MongoClient('mongodb://localhost:27017',
            username=db_obj["mongodbuser"],
            password=db_obj["mongodbpassword"],
            authSource=db_obj["mongodbname"],
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        dbh = client[db_obj["mongodbname"]]
        
        query_obj = {}
        if in_json["queryvalue"] != "":
            query = str(in_json["queryvalue"])
            id_query = str(in_json["queryvalue"])
            id_query = id_query.replace("GLYDS", "")
            id_query = id_query.replace("DSBCO_", "")
            cond_objs = []
            cond_objs.append({"bco_id":{'$regex':id_query, '$options': 'i'}})
            cond_objs.append({"provenance_domain.name":{'$regex':query, '$options': 'i'}})
            cond_objs.append({"provenance_domain.contributors.name":{'$regex':query, '$options': 'i'}})
            cond_objs.append({"io_domain.output_subdomain.uri.filename":{'$regex':query, '$options': 'i'}})
            query_obj = { "$or": cond_objs }

        #print query_obj


        out_json["datasets"] = []
        bco_collection = "c_bco_v-" + config_json["datarelease"]

        category_list = []
        for doc in dbh[bco_collection].find(query_obj).sort("bco_id", 1):
            if "extension_domain" not in doc:
                continue

            retired_flag = ""
            if "dataset_categories" in doc["extension_domain"]:
                for o in doc["extension_domain"]["dataset_categories"]:
                    cat_name = o["category_name"].replace(" ", "_")
                    cat_value = o["category_value"]
                    if cat_name in ["dataset_status", "status"] and cat_value == "retired":
                        retired_flag = "retired"
            if retired_flag == "retired":
                continue

            if "dataset_categories" in doc["extension_domain"]:
                for o in doc["extension_domain"]["dataset_categories"]:
                    cat_name = o["category_name"].replace(" ", "_")
                    if cat_name not in category_list and cat_name not in ["tags"]:
                        category_list.append(cat_name)
        out_json["categorylist"] = category_list

      

        doc_list = []
        for doc in dbh[bco_collection].find(query_obj).sort("bco_id", 1):
            if "bco_id" in doc:
                bco_id = doc["bco_id"].split("/")[-1]
                doc_list.append(doc)

        for doc in doc_list:
            doc.pop("_id")
            file_name = ""
            if "io_domain" in doc:
                if "output_subdomain" in doc["io_domain"]:
                    if doc["io_domain"]["output_subdomain"] != []:
                        file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"].strip()
            if file_name == "":
                continue
            
            category_dict = {}
            #for cat_name in category_list:
            #    category_dict[cat_name] = "Other"
            
            if "dataset_categories" in doc["extension_domain"]:
                for o in doc["extension_domain"]["dataset_categories"]:
                    cat_name = o["category_name"].replace(" ", "_")
                    cat_value = o["category_value"].strip()
                    if cat_value.strip() == "":
                        continue
                    category_dict[cat_name] = cat_value
                    #if cat_name == "species":
                    #    category_dict[cat_name] = cat_value[0].upper()
                    #    category_dict[cat_name] += cat_value.split(" ")[0][1:].lower() + " "
                    #    if len(cat_value.split(" ")) > 1:
                    #        category_dict[cat_name] += cat_value.split(" ")[1].lower()
            retired_flag = ""
            for s in ["dataset_status", "status"]:
                if s in category_dict:
                    if category_dict[s].lower() in ["retired"]:
                        retired_flag = "retired"
            if retired_flag == "retired":
                continue
            
            file_type = file_name.split(".")[-1]
            desc = doc["usability_domain"][0] if doc["usability_domain"] != [] else ""
            desc = desc[0:100] if len(desc) > 100 else desc
            obj = {
                "_id": doc["bco_id"].split("/")[-1],
                "title":doc["provenance_domain"]["name"],
                "description":desc,
                "filename": file_name,
                "filetype": file_type,
                "categories":category_dict
            }
            if file_type in ["csv", "tsv", "txt"]:
                header_row, body_rows = [], []
                delim = "," if file_type == "csv" else "\t"
                file_path = data_path + "/%s/reviewed/%s" %(current_release, file_name)
                if os.path.isfile(file_path):
                    with open(file_path, 'r') as FR:
                        data_frame = csv.reader(FR, delimiter=delim, quotechar='"')
                        row_count = 0
                        for row in data_frame:
                            row_count += 1
                            if len(row) < 1:
                                continue
                            if row_count == 1:
                                if len(row) == 1:
                                    header_row = [row[0][0:12] + "..."]
                                else:
                                    header_row = [row[0][0:12] + "...", row[1][0:12] + "..."]
                            elif row_count < 4:
                                if len(row) == 1:
                                    body_rows.append([row[0]])
                                else:
                                    body_rows.append([row[0], row[1][0:20] + "..."])
                            else:
                                break
                obj["minitable"] = {"content":body_rows, "headers":header_row, "colwidth": ["40%", "20%"]}
            else:
                #if "species" in category_dict:
                #    parts = category_dict["species"].split(" ")
                #    p1, p2 = parts[0], ""
                #    sp = p1[0].upper() + p1[1:]
                #    if len(parts) > 1:
                #        p2 = parts[1]
                #        sp += p2[0].upper() + p2[1:]
                #    obj["iconfilename"] = "%s_icon_%s.png" % (file_type, sp)
                obj["iconfilename"] = "%s_icon_all.png" % (file_type)
                if file_name.find("glycan_images") != -1:
                    image_type = file_name.split("_")[2]
                    obj["iconfilename"] = "%s_icon_all.png" % (image_type)
            out_json["datasets"].append(obj)

    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

