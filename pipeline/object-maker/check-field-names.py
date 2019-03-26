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

        field_dict = {}
        in_file =  "unreviewed/field_names.csv"
        libgly.load_sheet(field_dict, in_file, ",")

        field_list = []
        for row in field_dict["data"]:
            field = row[1] if row[1] != "" else row[0]
            field_list.append(field)


        pattern = "temp/*.csv"
        for in_file in glob.glob(pattern):
            file_name = in_file.split("/")[-1]
            if file_name == "field_names.csv":
                continue
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            for field in data_frame["fields"]:
                if field not in field_list:
                    print "undefined", field, file_name



if __name__ == '__main__':
	main()

