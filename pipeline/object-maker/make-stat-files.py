import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
from Bio import SeqIO

sys.path.append('../../glytools/')
import libgly


__version__="1.0"
__status__ = "Dev"



###############################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        data_dir = "reviewed/"

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


            in_file = "generated/misc/field_names.json"
            field_dict = json.loads(open(in_file,  "r").read())
            seen = {}
            seen_file = {}
            doc_list = list(dbh["c_bco"].find({}))
            for doc in doc_list:
                if "bco_id" not in doc:
                    continue
                if "extension_domain" not in doc:
                    continue
                bco_id = doc["bco_id"]
                status_list = []

                if "dataset_categories" in doc["extension_domain"]:
                    for o in doc["extension_domain"]["dataset_categories"]:
                        if o["category_name"] == "status":
                            status_list.append(o["category_value"].lower())
                if "retired" in status_list:
                    continue
                if "io_domain" not in doc:
                    continue
                if "output_subdomain" not in doc["io_domain"]:
                    continue
                if len(doc["io_domain"]["output_subdomain"]) < 1:
                    continue
                if "uri" not in doc["io_domain"]["output_subdomain"][0]:
                    continue
                if "filename" not in doc["io_domain"]["output_subdomain"][0]["uri"]:
                    continue
                file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"]
                if file_name == "":
                    print "Empty io_domain for bco %s!" % (bco_id)
                    continue
                file_name = file_name.strip()
                file_ext = file_name.split(".")[-1]
                in_file = data_dir + "/" + file_name
                if os.path.isfile(in_file) == True:
                    seen_file[file_name] = True
            
            for doc in dbh["c_bco"].find({}):
                if "bco_id" not in doc:
                    continue
                bco_id = doc["bco_id"]
                if "io_domain" not in doc:
                    continue
                if doc["io_domain"]["output_subdomain"] == []:
                    continue
               
                status_list = []
                if "dataset_categories" in doc["extension_domain"]:
                    for o in doc["extension_domain"]["dataset_categories"]:
                        if o["category_name"] == "status":
                            status_list.append(o["category_value"].lower())

                if "retired" in status_list:
                    continue

                bco_idx = doc["bco_id"].split("/")[-1]
                file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"].strip()
                if file_name not in seen_file:
                    print "** %s skipped because %s does not exist!" % (bco_idx,file_name)
                    continue

                
                stat_file_name = ".".join(file_name.split(".")[0:-1]) + ".stat.csv"
                stat_file = data_dir + "/" + stat_file_name
                flag = False
                for o in doc["io_domain"]["output_subdomain"]:
                    url = "http://data.glygen.org/ln2wwwdata/reviewed/" + file_name
                    o["uri"]["uri"] = url
                    if o["uri"]["filename"].split(".")[-2] == "stat":
                        flag = True


                #update io_domain
                if flag == False:
                    query_obj = {"bco_id":doc["bco_id"]}
                    update_obj = {"io_domain": doc["io_domain"]}
                    url = "http://data.glygen.org/ln2wwwdata/reviewed/" + stat_file_name
                    o = {
                        "mediatype": "csv", 
                        "uri":{"access_time":"", "sha1_chksum":"", "uri":url, "filename":stat_file_name}
                    }
                    update_obj["io_domain"]["output_subdomain"].append(o)
                    result = dbh["c_bco"].update_one(query_obj, {'$set': update_obj}, upsert=True)
                    print "added %s to %s" % (stat_file_name, doc["bco_id"])


                if file_name not in seen:
                    seen[file_name] = True
                    print "%s, %s " % (doc["bco_id"].split("/")[-1],stat_file.split("/")[-1])
                    FW = open(stat_file, "w")
                    row = ["unique_values","field_name", "field_description"]
                    FW.write("\"%s\"\n" % ("\",\"".join(row)))
                    if file_ext == "csv":
                        data_frame = {}
                        libgly.load_sheet(data_frame, in_file, ",")
                        stat_dict = {}
                        for row in data_frame["data"]:
                            for j in xrange(0, len(row)):
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
                            FW.write("\"%s\"\n" % ("\",\"".join(row)))
                    elif file_ext == "fasta":
                        n = len(list(SeqIO.parse(in_file, "fasta")))
                        row = [str(n), "sequence","Amino acid sequence"]
                        FW.write("\"%s\"\n" % ("\",\"".join(row)))
                    elif file_ext == "nt":
                        n = 0
                        with open(in_file, "r") as FR:
                            for line in FR:
                                n += 1
                        row = [str(n), "triple","RDF triple"]
                        FW.write("\"%s\"\n" % ("\",\"".join(row)))
                    FW.close()
                
        except pymongo.errors.ServerSelectionTimeoutError as err:
            return {}, {"error_list":[{"error_code": "open-connection-failed"}]}
        except pymongo.errors.OperationFailure as err:
            return {}, {"error_list":[{"error_code": "mongodb-auth-failed"}]}



if __name__ == '__main__':
	main()

