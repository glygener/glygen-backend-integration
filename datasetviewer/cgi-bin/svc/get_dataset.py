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
import json
import util

import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"



##############################
def get_preview(doc, path_obj, obj_ver):


    out_json = {"status":1, "errormsg":"", "info":{}}
    out_json["info"]["title"] = doc["provenance_domain"]["name"] 
    out_json["info"]["description"] = doc["provenance_domain"]["name"]

    file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"]
    file_type = file_name.split(".")[-1]


    out_json["info"]["filename"] = file_name
    out_json["info"]["filetype"] = file_type
    out_json["info"]["objid"] = doc["bco_id"].split("/")[-1].split(".")[0]
    out_json["info"]["filestatus"] = False               
    out_json["txtbuffer"] = ""

    
    cmd = "cat " + path_obj["htmlpath"] + "/ln2wwwdata/reviewed/release-notes.txt"
    release_info = commands.getoutput(cmd).strip()
    current_version = release_info.split(" ")[0]
    out_json["versions"] = [release_info]
    out_json["selectedversion"] = obj_ver if obj_ver != None else current_version

    file_list = glob.glob(path_obj["htmlpath"] + "/ln2wwwdata/reviewed/v-*/release-notes.txt")
    for rel_file in file_list:
        if os.path.exists(rel_file):
            cmd = "cat " + rel_file
            release_info = commands.getoutput(cmd).strip()
            ver = release_info.split(" ")[0]
            file_path = path_obj["htmlpath"] + "/ln2wwwdata/reviewed/%s/%s" % (ver,out_json["info"]["filename"])
            if os.path.isfile(file_path) == True:
                out_json["versions"].append(release_info)


    data_dir = path_obj["htmlpath"] + "/ln2wwwdata/reviewed/"
    if obj_ver != None and obj_ver != current_version:
            data_dir += obj_ver + "/"
    file_path = data_dir + out_json["info"]["filename"]
    if os.path.exists(file_path) == False:
        out_json["info"]["filestatus"] = False
        return out_json
    else:
        species_short = file_name.split("_")[0]
        if file_type == "csv":
            out_json["dataframe"] = []
            with open(file_path, 'r') as FR:
                csvGrid = csv.reader(FR, delimiter=',', quotechar='"')
                rowCount = 0
                for row in csvGrid:
                    rowCount += 1
                    out_json["dataframe"].append(row)
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
            in_file = path_obj["htmlpath"] + "/ln2wwwdata/reviewed/" + file_name
            for record in SeqIO.parse(in_file, "fasta"):
                seqCount += 1
                seqId = record.id
                seqDesc = record.description
                seqBody = str(record.seq.upper())
                out_json["seqobjects"].append({"seqid":seqId, "seqdesc":seqDesc, "seqbody":seqBody})
                if seqCount == 10:
                    break
        elif file_type in ["rdf"]:
            in_file = path_obj["htmlpath"] + "/ln2wwwdata/reviewed/" + file_name
            with open(in_file, 'r') as FR:
                lineCount = 0
                for line in FR:
                    lineCount += 1
                    out_json["txtbuffer"] += line
                    if lineCount == 100:
                        break
        elif file_type in ["gp", "gb", "nt"]:
            in_file = path_obj["htmlpath"] + "/ln2wwwdata/reviewed/" + file_name
            with open(in_file, 'r') as FR:
                lineCount = 0
                for line in FR:
                    lineCount += 1
                    out_json["txtbuffer"] += line
                    if lineCount == 200:
                        break
        elif file_type in ["png"]:
            file_list = glob.glob(path_obj["htmlpath"] + "/ln2wwwdata/glycanimages/*.png")
            for f in file_list[1000:1010]:
                recordId = f.split("/")[-1].split(".")[0]
                url = root_obj["htmlroot"] + "/ln2wwwdata/glycanimages/"  + f.split("/")[-1]
                out_json["txtbuffer"] +=  "%s <br>" % (recordId)
                out_json["txtbuffer"] += "<img src=\"%s\"><hr><br>" % (url)
        elif file_type in ["aln"]:
            file_list = glob.glob(path_obj["htmlpath"] + "/ln2wwwdata/aln/"+species_short+"/*.aln")
            for f in file_list[0:10]:
                recordId = f.split("/")[-1].split(".")[0]
                out_json["txtbuffer"] +=  "%s" % (recordId)
                with open(f, 'r') as FR:
                    for line in FR:
                        out_json["txtbuffer"] += line
                    out_json["txtbuffer"] += "\n\n"
        else:
            out_json["txtbuffer"] +=  "Please implement service for " + file_type + " preview!"


    return out_json



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-o","--objid",action="store",dest="objid",help="Object ID")
    parser.add_option("-v","--objver",action="store",dest="objver",help="Object Version")


    (options,args) = parser.parse_args()
    for key in ([options.objid]):
        if not (key):
            parser.print_help()
            sys.exit(0)


    global config_json
    global db_obj
    global client


    out_json = {}

    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        db_obj = config_json[config_json["server"]]["dbinfo"]

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
        options.objid = options.objid.replace("GLYDS", "DSBCO_")
        obj_id = "http://data.glygen.org/%s" % (options.objid)
        obj_ver = options.objver
        query_obj = {"bco_id":obj_id}
        bco_collection = "c_bco_v-" + config_json["datarelease"]
        doc = dbh[bco_collection].find_one(query_obj)
        out_json = get_preview(doc, config_json[config_json["server"]]["pathinfo"], obj_ver)


    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

