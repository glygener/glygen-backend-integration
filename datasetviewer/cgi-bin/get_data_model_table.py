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


__version__="1.0"
__status__ = "Dev"



#######################################
def main():

	usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)

        global PHASH
        global AUTH
        PHASH = {}

        custom_config_json = json.loads(open("conf/config.custom.json", "r").read())
        path_obj = custom_config_json["pathinfo"]

	outJson = { "taskstatus":1}
	datasetFile = path_obj["htmlpath"] + "/content/datamodel.csv"
        outJson["dataframe"] = []
	with open(datasetFile, 'r') as FR:
            csvGrid = csv.reader(FR, delimiter=',', quotechar='"')
            rowCount = 0
            for row in csvGrid:
                rowCount += 1
                #if rowCount > 6:
                #    continue
                outJson["dataframe"].append(row)
               
        #print outJson
        print json.dumps(outJson, indent=4)



if __name__ == '__main__':
        main()



