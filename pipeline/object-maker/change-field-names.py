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

        field_dict = {}
        in_file = data_dir + "/field_names.csv"
        libgly.load_sheet_as_dict(field_dict, in_file, ",", "field_name_current")
        
        pattern = data_dir + "/*.csv"
        for in_file in glob.glob(pattern):
            file_name = in_file.split("/")[-1]
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            for j in xrange(0, len(data_frame["fields"])):
                field = data_frame["fields"][j]
                if field in field_dict["data"]:
                    if field_dict["data"][field][0][0] != "":
                        data_frame["fields"][j] = field_dict["data"][field][0][0]
            
            out_file = "temp/%s" % (file_name)
            FW = open(out_file, "w")
            FW.write("\"%s\"\n" % ( "\",\"".join(data_frame["fields"])))
            for row in data_frame["data"]:
                FW.write("\"%s\"\n" % ( "\",\"".join(row)))
            FW.close()



if __name__ == '__main__':
	main()

