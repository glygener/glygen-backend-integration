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




###############################
def main():


    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    msg = "Input JSON text"
    parser.add_option("-j","--in_json",action="store",dest="in_json",help=msg)

    (options,args) = parser.parse_args()
    for file in ([options.in_json]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    global config_json
    global db_obj
    global path_obj
    global client


    out_json = {}
    in_json  = json.loads(options.in_json)

    species_dict = {"human":"Homo sapiens", "mouse":"Mus musculus"}
    cat_dict = {"glycan":"Glycan centric", "protein":"Protein centric", "proteoform":"Proteoform centric"}
  

    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        db_obj = config_json[config_json["server"]]["dbinfo"]
        path_obj = config_json[config_json["server"]]["pathinfo"]
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
        
        query_obj = {}
        if in_json["queryvalue"] != "":
            query = str(in_json["queryvalue"])
            cond_objs = []
            cond_objs.append({"bco_id":{'$regex':query, '$options': 'i'}})
            cond_objs.append({"provenance_domain.name":{'$regex':query, '$options': 'i'}})
            cond_objs.append({"provenance_domain.contributors.name":{'$regex':query, '$options': 'i'}})
            query_obj = { "$or": cond_objs }
        
        out_json["datasets"] = []
        bco_collection = "c_bco_v-" + config_json["datarelease"]
        for doc in dbh[bco_collection].find(query_obj).sort("bco_id", 1):
            doc.pop("_id")
            file_name = ""
            if "io_domain" in doc:
                if "output_subdomain" in doc["io_domain"]:
                    if doc["io_domain"]["output_subdomain"] != []:
                        file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"]
            if file_name == "":
                continue

            file_type = file_name.split(".")[-1]
            species_id = file_name.split("_")[0]
            species = species_dict[species_id] if species_id in species_dict else "Other"
            category_id = file_name.split("_")[1]
            category = cat_dict[category_id] if category_id in cat_dict else "Other"
            obj = {
                "_id": doc["bco_id"].split("/")[-1],
                "title":doc["provenance_domain"]["name"],
                "description":doc["provenance_domain"]["name"],
                "filename": file_name,
                "filetype": file_type,
                "status": 1,
                "category": category,
                "species": species
            }
            if file_type == "csv":
                header_row, body_rows = [], []
                file_path = path_obj["htmlpath"] + "ln2wwwdata/reviewed/" +  file_name
                if os.path.isfile(file_path):
                    with open(file_path, 'r') as FR:
                        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
                        row_count = 0
                        for row in data_frame:
                            row_count += 1
                            if len(row) < 2:
                                continue
                            if row_count == 1:
                                header_row = [row[0][0:12] + "...", row[1][0:12] + "..."]
                            elif row_count < 4:
                                body_rows.append([row[0], row[1][0:20] + "..."])
                            else:
                                break
                obj["minitable"] = {"content":body_rows, "headers":header_row, "colwidth": ["40%", "20%"]}

            else:
                obj["iconfilename"] = "%s_icon_%s.png" % (file_type, species_id)
            out_json["datasets"].append(obj)
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

