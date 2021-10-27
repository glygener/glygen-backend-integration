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
import gzip

import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"

def get_top_lines(in_file):

    buffer = ""
    FR = gzip.open(in_file, 'rb') if in_file[-3:] == ".gz" else open(in_file, 'r')
    l_count = 0 
    for line in FR: 
        l_count += 1
        buffer += line
        if l_count == 1000:
            break
    FR.close()
    return buffer


def check_bco_existence(bco_id, rel):
    
    query_obj = {"bco_id":bco_id}
    bco_collection = "c_bco_v-" + rel
    doc = dbh[bco_collection].find_one(query_obj)
    return  doc != None


##############################
def get_preview(doc, path_obj, obj_ver):



    out_json = {"status":1, "errormsg":"", "info":{}}
    out_json["info"]["title"] = doc["provenance_domain"]["name"] 
    desc = doc["usability_domain"][0] if doc["usability_domain"] != [] else ""
    #desc = desc[0:100] if len(desc) > 100 else desc
    out_json["info"]["description"] =  desc


    file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"].strip()
    file_type = file_name.split(".")[-1]


    out_json["info"]["filename"] = file_name
    out_json["info"]["filetype"] = file_type
    out_json["info"]["rndrtype"] = "html"
    out_json["info"]["objid"] = doc["bco_id"].split("/")[-1].split(".")[0]
    out_json["info"]["filestatus"] = False               
    out_json["txtbuffer"] = ""


    rel_data_path = data_path + "/%s/" % (selected_release)
    rel_data_root = data_root + "/%s/" % (selected_release)

    rel_file = data_path + "/%s/releaseinfo/all_history.json" % (current_release)

    rel_json = json.loads(open(rel_file, "r").read())
    id_list = rel_json.keys()
    rel_date_dict = {}
    for bco_id in id_list:
        new_bco_id = bco_id.replace("DSBCO_", "GLY_")
        rel_json[new_bco_id] = rel_json[bco_id]
        for rel in rel_json[bco_id]:
            rel_date_dict[rel] = rel_json[bco_id][rel]["releasedate"]
        
    bco_id = doc["bco_id"].split("/")[-1]
    out_json["versions"] = []
    out_json["selectedversion"] = selected_release
    r_list = [current_release.split("-")[-1]]
    if bco_id in rel_json:
        r_list = sorted(rel_json[bco_id], reverse=True)

    for rel in r_list:
        rel_date = rel_date_dict[rel]
        bco_exists = check_bco_existence(doc["bco_id"], rel)
        if bco_exists == True:
            out_json["versions"].append("v-" + rel + " " + rel_date)


    
    file_path = rel_data_path + "reviewed/"+  out_json["info"]["filename"]
    if os.path.exists(file_path) == False:
        out_json["info"]["filestatus"] = False
        return out_json
    else:
        species_short = file_name.split("_")[0]
        if file_type in ["csv", "tsv"]:
            delim = "," if file_type == "csv" else "\t"
            out_json["dataframe"] = []
            with open(file_path, 'r') as FR:
                csvGrid = csv.reader(FR, delimiter=delim, quotechar='"')
                rowCount = 0
                for row in csvGrid:
                    rowCount += 1
                    tmp_row = []
                    for val in row:
                        val = str(val).encode('utf-8').strip()
                        tmp_row.append(val)
                    out_json["dataframe"].append(tmp_row)
                    if rowCount == 1:
                        tmpList = []
                        for j in xrange(0, len(row)):
                            tmpList.append("string")
                        out_json["dataframe"].append(tmpList)
                    if rowCount == 85:
                        break
        elif file_type == "fasta":
            out_json["seqobjects"] = []
            seqCount = 0
            in_file = rel_data_path + "/reviewed/" + file_name
            for record in SeqIO.parse(in_file, "fasta"):
                seqCount += 1
                seqId = record.id
                seqDesc = record.description
                seqBody = str(record.seq.upper())
                out_json["seqobjects"].append({"seqid":seqId, "seqdesc":seqDesc, "seqbody":seqBody})
                if seqCount == 10:
                    break
        elif file_type in ["gz"] and file_name.find("glycan_images") != -1:
            i_type = file_name.split("_")[2]
            i_path = rel_data_path + "/glycanimages_%s/" % (i_type)
            i_root = rel_data_root + "/glycanimages_%s/" % (i_type)
            file_list = glob.glob(i_path + "*.png")
            for f in file_list[1000:1010]:
                recordId = f.split("/")[-1].split(".")[0]
                url = i_root +  f.split("/")[-1]
                out_json["txtbuffer"] +=  "%s <br>" % (recordId)
                out_json["txtbuffer"] += "<img src=\"%s\"><hr><br>" % (url)
        elif file_type.lower() in ["gif", "jpeg", "jpg"]:
            url = rel_data_root + "/reviewed/%s" % (file_name)
            out_json["txtbuffer"] += "<img src=\"%s\"><br>" % (url)
        elif file_type.lower() in ["mp4"]:
            url = rel_data_root + "/reviewed/%s" % (file_name)
            out_json["txtbuffer"] += "<video controls=\"controls\"><source src=\"%s\" type=\"video/mp4\"></video><br>" % (url)
        elif file_type in ["gz","gp", "gb", "nt"]:
            out_json["info"]["rndrtype"] = "text"
            in_file = rel_data_path + "/reviewed/" + file_name
            out_json["txtbuffer"] += get_top_lines(in_file)
        else:
            out_json["txtbuffer"] +=  "Please implement service for " + file_type + " preview!"


    return out_json



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



    global dbh
    global config_json
    global db_obj
    global client
    global root_obj
    global data_path
    global data_root
    global selected_release
    global current_release


    print "Content-Type: application/json"
    print   

    

    out_json = {}
    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        custom_config_json = json.loads(open("conf/config.custom.json", "r").read())
        db_obj = custom_config_json[config_json["server"]]["dbinfo"]
        path_obj = custom_config_json[config_json["server"]]["pathinfo"]
        root_obj = custom_config_json[config_json["server"]]["rootinfo"]
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
        
        obj_id = in_json["objid"].replace(db_obj["dsprefix"], "")
        obj_id = db_obj["bcourl"] % (obj_id)
        obj_ver = in_json["objver"] if "objver" in in_json else ""
        query_obj = {"bco_id":obj_id}


        current_release = "v-" + config_json["datarelease"]
        selected_release = obj_ver if obj_ver not in [None, ""] else current_release
        data_path = path_obj["htmlpath"] + "/ln2releases/"
        data_root = root_obj["htmlroot"] + "/ln2releases/"

        bco_collection = "c_bco_" + selected_release
        doc = dbh[bco_collection].find_one(query_obj)

        if doc != None:
            out_json = get_preview(doc, custom_config_json[config_json["server"]]["pathinfo"], obj_ver)

    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    version_one, version_two = config_json["moduleversion"],config_json["datarelease"]
    out_json["dmversions"] = "Viewer-%s | Data-%s" % (version_one, version_two)
    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
        main()

