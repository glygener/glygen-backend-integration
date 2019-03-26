import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient

sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"



###############################
def main():

        config_obj = json.loads(open("conf/config.json", "r").read())
        db_obj = config_obj[config_obj["server"]]["dbinfo"]
        data_dir = "unreviewed-backup/"

        sheet_obj = {}
        in_file = data_dir + "/field_names.csv"
        libgly.load_sheet(sheet_obj, in_file, ",")
        
        correction_dict = {}
        for row in sheet_obj["data"]:
            if row[1] != "" and row[0] != row[1]:
                correction_dict[row[0]] = row[1]


        pattern = "*.py"
        for in_file in glob.glob(pattern):
            if in_file == "check-field-names-in-code.py":
                continue
            with open(in_file, "r") as FR:
                for line in FR:
                    if "index(" in line:
                        f = line.split("index(")[1].split(")")[0].replace("\"", "")
                        if f in correction_dict:
                            print in_file, f, correction_dict[f]




if __name__ == '__main__':
	main()

