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
		outputFile = serverJson["pathinfo"]["htmlpath"] + "/jsondb/commentsdb.json"
		commentsJson = json.loads(open(outputFile, "r").read())
		if inJson["objid"] in commentsJson:
			outJson["comments"] = commentsJson[inJson["objid"]]
	except Exception, e:
                outJson["taskStatus"] = 0
                outJson["errorMsg"] = errorMsg if errorMsg != "" else str(e)
	
	print json.dumps(outJson, indent=4, separators=(',',':'))


if __name__ == '__main__':
        main()



