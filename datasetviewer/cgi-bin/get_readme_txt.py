#!/usr/bin/python
import os,sys
import string
import commands
import csv
import traceback

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import cgi
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
    global client
    global current_release
    global data_path
    global data_root

    print "Content-Type: text/html"
    print
    
    #out_buffer += json.dumps(in_json, indent=4)
    #sys.exit()


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
        out_json = {"taskstatus":0, "errormsg":"Loading config failed!"}

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
        obj_id = in_json["objid"].replace(db_obj["dsprefix"], "")
        obj_id = db_obj["bcourl"] % (obj_id)
        obj_ver = in_json["objver"] if "objver" in in_json else ""
        query_obj = {"bco_id":obj_id}
        bco_collection = "c_bco_v-" + config_json["datarelease"]
        doc = dbh[bco_collection].find_one(query_obj)

        if doc == None:
            out_buffer += "\n\n\tReadme file does not exist for object %s!" % (obj_id)
        else:
            creator_list = []
            author_list = []
            for o in doc["provenance_domain"]["contributors"]:
                for p in ["createdBy"]:
                    if p in o["contribution"]:
                        creator_list.append("%s [%s]" % (o["name"], o["email"]))
                for p in ["contributedBy", "authoredBy"]:
                    if p in o["contribution"]:
                        author_list.append("%s [%s]" % (o["name"], o["email"]))
            
            modified = ""
            if "modified" in doc["provenance_domain"]:
                modified = doc["provenance_domain"]["modified"]
            checksum = doc["checksum"] if "checksum" in doc else ""

            out_buffer = ""
            out_buffer += "<b>Name</b>\n   - %s" % (doc["provenance_domain"]["name"])
            out_buffer += "\n\n<b>Version</b>\n   - %s" % (doc["provenance_domain"]["version"])
            out_buffer += "\n\n<b>Created on</b>\n  - %s" % (doc["provenance_domain"]["created"])
            out_buffer += "\n\n<b>Modified on</b>\n   - %s" % (modified)
            out_buffer += "\n\n<b>Created by</b>\n   - %s" % (", ".join(creator_list))
            out_buffer += "\n\n<b>Authors</b>\n   - %s" % (", ".join(author_list))
            out_buffer += "\n\n<b>Checksum</b>\n   - %s" % (checksum)
            if doc["usability_domain"]  != []:
                out_buffer += "\n\n<b>Usability Domain</b>\n   - %s" % ("\n   - ".join(doc["usability_domain"]))
            
            if doc["description_domain"]["keywords"] != [] or doc["description_domain"]["pipeline_steps"] != []:
                out_buffer += "\n\n<b>Description Domain</b>"
                if "keywords" in doc["description_domain"]:
                    if doc["description_domain"]["keywords"] != []:
                        out_buffer += "\n   Keywords\n   - %s" % ("\n   - ".join(doc["description_domain"]["keywords"]))
                if "pipeline_steps" in doc["description_domain"]:
                    if doc["description_domain"]["pipeline_steps"] != []:
                        out_buffer += "\n\n   Pipeline Steps"
                        for o in doc["description_domain"]["pipeline_steps"]:
                            out_buffer += "\n   - Step-%s: %s" % (o["step_number"], o["description"])


            if doc["execution_domain"] != {}:
                out_buffer += "\n\n<b>Execution Domain</b>"
                if "script" in doc["execution_domain"]:
                    out_buffer += "\n   Scripts"
                    for o in doc["execution_domain"]["script"]:
                        out_buffer += "   - %s\n" % (o["uri"]["uri"])
                if "software_prerequisites" in doc["execution_domain"]:
                    out_buffer += "\n\n   Software Prerequisites"
                    for o in doc["execution_domain"]["software_prerequisites"]:
                        out_buffer += "\n      Name: %s" % (o["name"])
                        out_buffer += "\n      Version: %s" % (o["version"])
                        out_buffer += "\n      URI: %s" % (o["uri"]["uri"])
                        out_buffer += ""

            if doc["io_domain"] != {}:
                out_buffer += "\n\n<b>I/O Domain</b>"
                if "input_subdomain" in doc["io_domain"]:
                    if doc["io_domain"]["input_subdomain"] != []:
                        out_buffer += "\n   Input Subdomain"
                        for o in doc["io_domain"]["input_subdomain"]:
                            out_buffer += "\n      Name: %s" % (o["uri"]["filename"])
                            #out_buffer += "      Media Type: %s" % (o["mediatype"])
                            out_buffer += "\n      URI: %s" % (o["uri"]["uri"])
                            out_buffer += ""
                if "output_subdomain" in doc["io_domain"]:
                    if doc["io_domain"]["output_subdomain"] != []:
                        out_buffer += "\n   Output Subdomain"
                        for o in doc["io_domain"]["output_subdomain"]:
                            out_buffer += "\n      Name: %s" % (o["uri"]["filename"])
                            #out_buffer += "      Media Type: %s" % (o["mediatype"])
                            out_buffer += "\n      URI: %s" % (o["uri"]["uri"])
                            out_buffer += ""
                            if o["uri"]["filename"].split(".")[-2] == "stat":
                                stat_file = data_path +"/%s/reviewed/%s"
                                stat_file = stat_file % (current_release, o["uri"]["filename"])
                                if os.path.isfile(stat_file):
                                    stat_lines = open(stat_file, "r").read().split("\n")
                                    out_buffer += "\n\n<b>Fields</b>"
                                    for stat_line in stat_lines[1:]:
                                        if stat_line.strip() == "":
                                            continue
                                        s_p = stat_line.split(",")
                                        out_buffer += "\n      Name: %s" % (s_p[1])
                                        out_buffer += "\n     Description: %s" % (s_p[2])

                                        out_buffer += "\n      Unique Values: %s" % (s_p[0])
                                        out_buffer += ""   
        print "<pre>" 
        print out_buffer.encode(encoding='UTF-8',errors='strict')
        print "</pre>"
        sys.exit()
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print  json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

