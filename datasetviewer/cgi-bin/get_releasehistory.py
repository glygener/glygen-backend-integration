#!/usr/bin/python
import os,sys
reload(sys)
sys.setdefaultencoding("ISO-8859-1")

import string
import commands
import csv
import traceback

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import util
import cgi


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
    global client
    global root_obj
    global current_release
    global data_path
    global data_root

    print "Content-Type: application/json"
    print   

    

    out_json = {}
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
 
        
        rel_file = data_path + "/%s/releaseinfo/all_history.json" % (current_release)
        rel_json = json.loads(open(rel_file, "r").read())


        seen = {}
        bcoid2filenames = {}
        for bco_id in rel_json:
            for ver in rel_json[bco_id]:
                seen[ver] = True
                if bco_id not in bcoid2filenames:
                    bcoid2filenames[bco_id] = []
                if ver in rel_json[bco_id][ver]:
                    file_name = rel_json[bco_id][ver][ver]["filename"]
                    if file_name not in bcoid2filenames[bco_id]:
                        bcoid2filenames[bco_id].append(file_name)
        ver_list = sorted(seen.keys())

        cat_dict = {}
        selected_bco_list = []
        for ver in ver_list:
            bco_collection = "c_bco_v-" + ver
            for doc in dbh[bco_collection].find({}):
                if "bco_id" not in doc:
                    continue
                if "io_domain" not in doc:
                    continue
                if "output_subdomain" not in doc["io_domain"]:
                    continue
                if doc["io_domain"]["output_subdomain"] == []:
                    continue
                bco_id = doc["bco_id"].split("/")[-1]
                if "dataset_categories" in doc["extension_domain"]:
                    for o in doc["extension_domain"]["dataset_categories"]:
                        cat_name = o["category_name"].strip().lower()
                        cat_value = o["category_value"].strip().lower()
                        if cat_name != "" and cat_value != "" and cat_name != "tags":
                            if cat_name not in cat_dict:
                                cat_dict[cat_name] = []
                            if cat_value not in cat_dict[cat_name]:
                                cat_dict[cat_name].append(cat_value)

                        if in_json != {} and cat_name == in_json["category_name"] and cat_value == in_json["category_value"]:
                            selected_bco_list.append(bco_id)



        row_one, row_two = ["File Name"], ["string"]
        for ver in ver_list:
            row_one.append("ver-"+ver)
            row_two.append("number")
        row_one.append("")
        row_two.append("string")
        data_frame = [row_one, row_two]
        for bco_id in rel_json:
            bco_id_new = bco_id.replace("DSBCO_", db_obj["bcoprefix"])
            if selected_bco_list != [] and bco_id not in selected_bco_list:
                continue
            row = [bco_id + " (%s)" %(", ".join(bcoid2filenames[bco_id]))]
            for ver in ver_list:
                val = 0
                if ver in rel_json[bco_id]:
                    if ver in rel_json[bco_id][ver]:
                        if "recordcount" in rel_json[bco_id][ver][ver]:
                            val = rel_json[bco_id][ver][ver]["recordcount"]
                row.append(val)
            detail_link = "<a href=\"%s/history\" target=_>details</a>" % (bco_id_new)
            row.append(detail_link)
            data_frame.append(row)
        out_json = {"dataframe":data_frame, "categories":cat_dict, "taskstatus":1, "errormsg":""}
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    version_one, version_two = config_json["moduleversion"],config_json["datarelease"]
    out_json["dmversions"] = "Viewer-%s | Data-%s" % (version_one, version_two)
    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

