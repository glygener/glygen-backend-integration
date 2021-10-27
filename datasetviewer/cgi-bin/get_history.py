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
    global path_obj
    global current_release
    global data_path
    global data_root

    print "Content-Type: text/html"
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
        obj_id = in_json["objid"]


        file_path = data_path + "/%s/releaseinfo/all_history.json" % (current_release)
        if os.path.isfile(file_path) == False:
            print "\n\n\tNo history found for %s!" % (obj_id)
        else:
            rel_json = json.loads(open(file_path, "r").read())    
            fname_list = []
            obj_id_old = obj_id.replace(db_obj["bcoprefix"], "DSBCO_")
            his_obj = rel_json[obj_id] if obj_id in rel_json else rel_json[obj_id_old]
            for ver in his_obj:
                file_name = his_obj[ver][ver]["filename"]
                if file_name not in fname_list:
                    fname_list.append(file_name)
            print "<b>BCO ID</b> : %s<br>" % (obj_id)
            print "<b>File Name(s)</b> : %s<br><hr>" % (", ".join(fname_list))
            for ver in sorted(his_obj, reverse=True):
                o = his_obj[ver][ver]
                id_list_one = sorted(o["additions"]) if "additions" in o else []
                id_list_two = sorted(o["deletions"]) if "deletions" in o else []
                n_one = len(id_list_one)
                n_two = len(id_list_two)
                f_str = "<br><br><b>Version-%s</b> (%s records, %s additions, %s deletions)<br>"
                print f_str % (ver, o["recordcount"], n_one, n_two)
                s = "padding:10px;width:100%;background:#eee;border:0px solid #ccc;"
                nrows = "10" if n_one > 0 else "2"
                l_one = "\n".join(id_list_one)
                box_one = "<textarea rows=10 style=\"%s\">%s</textarea>" % (s, l_one)
                l_two = "\n".join(id_list_two)
                box_two = "<textarea rows=10 style=\"%s\">%s</textarea>" % (s, l_two)
                print "<table border=0>"
                print "<tr><td><i>Additions</i></td><td><i>Deletions</i></td></tr>"
                print "<tr><td>%s</td><td>%s</td></tr>" % (box_one, box_two)
                print "</table>"
        sys.exit()
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

