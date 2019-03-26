#!/usr/bin/python
import os,sys
import string
import cgi
import commands
from optparse import OptionParser
import csv
import json
import util
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"



#######################################
def main():

	usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        msg = "Input JSON text"
        parser.add_option("-j","--in_json",action="store",dest="in_json",help=msg)

        global PHASH
        global AUTH
        PHASH = {}

        (options,args) = parser.parse_args()
        for file in ([options.in_json]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

        config_obj = json.loads(open("conf/config.json", "r").read())
        serverJson  =  config_obj[config_obj["server"]]
        path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
        root_obj =  config_obj[config_obj["server"]]["rootinfo"]
        db_obj = config_obj[config_obj["server"]]["dbinfo"]



        in_json  = json.loads(options.in_json)
	out_json = {"injson":in_json, "taskstatus":1}

        client = MongoClient('mongodb://localhost:27017')
        db = client[db_obj["dbname"]]
        coll = "c_metadata"
        obj_id = int(in_json["objectid"])

        query_obj = {"_id":obj_id}
        doc = db[coll].find_one(query_obj)

        species_short = "human" if doc["species"] ==  "Homo sapiens" else "mouse"
	
        if doc["filetype"] == "csv":
		out_json["dataframe"] = []
		inFile = serverJson["pathinfo"]["htmlpath"] + "/ln2wwwdata/reviewed/" + doc["filename"]
		with open(inFile, 'r') as FR:
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
                                if rowCount == 10:
					break
	elif doc["filetype"] == "fasta":
        	out_json["seqobjects"] = []
		seqCount = 0
		inFile = serverJson["pathinfo"]["htmlpath"] + "/ln2wwwdata/reviewed/" + doc["filename"]
		for record in SeqIO.parse(inFile, "fasta"):
                	seqCount += 1
			seqId = record.id
                	seqDesc = record.description
			seqBody = str(record.seq.upper())
			out_json["seqobjects"].append({"seqid":seqId, "seqdesc":seqDesc, "seqbody":seqBody})
			if seqCount == 10:
				break
	elif doc["filetype"] in ["rdf"]:
                out_json["rdflines"] = []
		inFile = serverJson["pathinfo"]["htmlpath"] + "/ln2wwwdata/reviewed/" + doc["filename"]
                with open(inFile, 'r') as FR:
                    lineCount = 0
                    for line in FR:
                        lineCount += 1
			print line.strip()
                        if lineCount == 100:
                            break
		sys.exit()
        elif doc["filetype"] in ["gp", "gb", "nt"]:
                inFile = serverJson["pathinfo"]["htmlpath"] + "/ln2wwwdata/reviewed/" + doc["filename"]
                with open(inFile, 'r') as FR:
                    lineCount = 0
                    for line in FR:
                        lineCount += 1
                        print line.replace("\n", "")
                        if lineCount == 200:
                            break
                sys.exit()
        elif doc["filetype"] in ["png"]:
            fileList = glob.glob(serverJson["pathinfo"]["htmlpath"] + "/ln2wwwdata/glycanimages/*.png")
            for f in fileList[0:10]:
                recordId = f.split("/")[-1].split(".")[0]
                url = serverJson["rootinfo"]["htmlroot"] + "/ln2wwwdata/glycanimages/"  + f.split("/")[-1]
                print "%s <br>" % (recordId)
                print "<img src=\"%s\"><hr><br>" % (url)
            sys.exit()
        elif doc["filetype"] in ["aln"]:
            fileList = glob.glob(serverJson["pathinfo"]["htmlpath"] + "/ln2wwwdata/aln/"+species_short+"/*.aln")
            for f in fileList[0:10]:
                recordId = f.split("/")[-1].split(".")[0]
                print "%s" % (recordId)
                with open(f, 'r') as FR:
                    for line in FR:
                        print line[:-1]
                    print "\n\n"
            sys.exit()
        else:
		out_json["errormsg"] = "Please implement service for " + doc["filetype"] + " preview!"
                

	print json.dumps(out_json, indent=4, separators=(',',':'))





if __name__ == '__main__':
        main()



