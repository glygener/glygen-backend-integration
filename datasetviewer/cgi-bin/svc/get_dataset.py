import os,sys
import string
import commands
import csv

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"




###############################
def main():


        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog version___")
        parser.add_option("-o","--objid",action="store",dest="objid",help="Object ID")
        
        (options,args) = parser.parse_args()
        for key in ([options.objid]):
            if not (key):
                parser.print_help()
                sys.exit(0)



        config_obj = json.loads(open("conf/config.json", "r").read())
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        
	client = MongoClient('mongodb://localhost:27017')
	db = client[db_obj["dbname"]]
        coll = "c_metadata"
        obj_id = int(options.objid.replace("GLYDS", ""))
       
        query_obj = {"_id":obj_id}
        doc = db[coll].find_one(query_obj)
        
        out_json = {"info":{}}
        out_json["info"]["title"] = doc["title"]
        out_json["info"]["description"] = doc["description"]
        out_json["info"]["filename"] = doc["filename"]
        out_json["info"]["filetype"] = doc["filetype"]

        out_json["info"]["objid"] = options.objid
        out_json["info"]["filestatus"] = False               
        out_json["txtbuffer"] = ""


        data_file = path_obj["htmlpath"] + "/datasets/reviewed/" + doc["filename"]
        if os.path.exists(data_file):
            out_json["info"]["filestatus"] = True

        species_short = "human" if doc["species"] ==  "Homo sapiens" else "mouse"
        
        if doc["filetype"] == "csv":
            out_json["dataframe"] = []
            with open(data_file, 'r') as FR:
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
                    if rowCount == 100:
                        break
        elif doc["filetype"] == "fasta":
                out_json["seqobjects"] = []
                seqCount = 0
                in_file = path_obj["htmlpath"] + "/datasets/reviewed/" + doc["filename"]
                for record in SeqIO.parse(in_file, "fasta"):
                        seqCount += 1
                        seqId = record.id
                        seqDesc = record.description
                        seqBody = str(record.seq.upper())
                        out_json["seqobjects"].append({"seqid":seqId, "seqdesc":seqDesc, "seqbody":seqBody})
                        if seqCount == 10:
                                break
        elif doc["filetype"] in ["rdf"]:
                in_file = path_obj["htmlpath"] + "/datasets/reviewed/" + doc["filename"]
                with open(in_file, 'r') as FR:
                    lineCount = 0
                    for line in FR:
                        lineCount += 1
                        out_json["txtbuffer"] += line
                        if lineCount == 100:
                            break
        elif doc["filetype"] in ["gp", "gb", "nt"]:
                in_file = path_obj["htmlpath"] + "/datasets/reviewed/" + doc["filename"]
                with open(in_file, 'r') as FR:
                    lineCount = 0
                    for line in FR:
                        lineCount += 1
                        out_json["txtbuffer"] += line
                        if lineCount == 200:
                            break
        elif doc["filetype"] in ["png"]:
            fileList = glob.glob(path_obj["htmlpath"] + "/datasets/glycanimages/*.png")
            for f in fileList[1000:1010]:
                recordId = f.split("/")[-1].split(".")[0]
                url = root_obj["htmlroot"] + "/datasets/glycanimages/"  + f.split("/")[-1]
                out_json["txtbuffer"] +=  "%s <br>" % (recordId)
                out_json["txtbuffer"] += "<img src=\"%s\"><hr><br>" % (url)
        elif doc["filetype"] in ["aln"]:
            fileList = glob.glob(path_obj["htmlpath"] + "/datasets/aln/"+species_short+"/*.aln")
            for f in fileList[0:10]:
                recordId = f.split("/")[-1].split(".")[0]
                out_json["txtbuffer"] +=  "%s" % (recordId)
                with open(f, 'r') as FR:
                    for line in FR:
                        out_json["txtbuffer"] += line
                    out_json["txtbuffer"] += "\n\n"
        else:
                out_json["txtbuffer"] +=  "Please implement service for " + doc["filetype"] + " preview!"
                

        print json.dumps(out_json, indent=4)
                        


if __name__ == '__main__':
	main()

