#!/usr/bin/python
import os,sys
import string
import cgi
import commands
from optparse import OptionParser
import json
import util
import datetime
import time




__version__="1.0"
__status__ = "Dev"



#~~~~~~~~~~~~~~~~~~~~~
def main():

   	usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        msg = "Input JSON text"
        parser.add_option("-j","--inJson",action="store",dest="inJson",help=msg)

	global PHASH
	global AUTH
	PHASH = {}


	(options,args) = parser.parse_args()
                
	for file in ([options.inJson]):
		if not (file):
			parser.print_help()
			sys.exit(0)
	

	configJson = json.loads(open("conf/config.json", "r").read())
	serverJson  =  configJson[configJson["server"]]
	
	inJson  = json.loads(options.inJson)
	outJson = {}

	errorMsg = ''
	try:
                      
		outJson = {"taskStatus":1, "inJson":inJson}

		#outputFile ="/tmp/comments.json"
		outputFile = serverJson["pathinfo"]["htmlpath"] + "/jsondb/commentsdb.json"
		commentsJson = json.loads(open(outputFile, "r").read())

		obj = {}		
		obj["createdts"] = datetime.datetime.now().strftime("%b %d, %Y %H:%M:%S")
		obj["fullname"] = inJson["fullname"]
		obj["email"] = inJson["email"]
		obj["comment"] = inJson["comment"]
		
		if inJson["objid"] not in commentsJson:
			commentsJson[inJson["objid"]] = []
		commentsJson[inJson["objid"]].append(obj);

		FW = open(outputFile, "w")
		FW.write("%s" % (json.dumps(commentsJson, indent=4)))
        	FW.close()
		cmd = "cp " + outputFile + " " + outputFile + "-backup"
       		x = commands.getoutput(cmd)
	except Exception, e:
                outJson["taskStatus"] = 0
                outJson["errorMsg"] = errorMsg if errorMsg != "" else str(e)
	
	print json.dumps(outJson, indent=4, separators=(',',':'))


if __name__ == '__main__':
        main()



