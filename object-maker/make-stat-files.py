import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
from Bio import SeqIO

import libgly


__version__="1.0"
__status__ = "Dev"


def get_stat_rows(in_file, field_dict):
    file_ext = in_file.split(".")[-1]
    row_list = []
    if file_ext in ["csv", "tsv"]:
        sep = "," if file_ext == "csv" else "\t"
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, sep)
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

    return row_list



def update_io_domain(doc):

    bco_id = doc["object_id"].split("/")[-2]
    in_obj_list = []
    for obj in doc["io_domain"]["output_subdomain"]:
        out_file_name = obj["uri"]["uri"].split("/")[-1]


    return



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
            for doc in dbh["c_bco"].find({}):
                stat_file_name = ".".join(file_name.split(".")[0:-1]) + ".stat.csv"
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

