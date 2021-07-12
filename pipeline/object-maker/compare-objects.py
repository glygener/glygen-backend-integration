import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-d","--dbtype",action="store",dest="dbtype",help="proteindb/glycandb")
    parser.add_option("-v","--ver",action="store",dest="ver",help="1.5.23")

    (options,args) = parser.parse_args()
    for key in ([options.dbtype, options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)


    dbtype = options.dbtype
    ver = "v-" + options.ver
    old_dir = "/data/projects/glygen/releases/data/%s/" % (ver)

    seen = {}
    file_list = glob.glob("jsondb/%s/*.json" % (dbtype))
    for file_two in file_list:
        canon = file_two.split("/")[-1].split(".")[0]
        in_file_one = old_dir + "jsondb/%s/%s.json" % (dbtype,canon)
        in_file_two = "jsondb/%s/%s.json" % (dbtype, canon)
        if os.path.isfile(in_file_one) == False or os.path.isfile(in_file_two) == False:
            continue
        doc_one = json.loads(open(in_file_one, "r").read())
        doc_two = json.loads(open(in_file_two, "r").read())
        for key in doc_one:
            str_one = json.dumps(doc_one[key])
            str_two = json.dumps(doc_two[key])
            if str_one != str_two:
                combo_id = key
                if dbtype == "proteindb":
                    tax_name = doc_one["species"][0]["name"]
                    combo_id = "%s %s" % (key, tax_name)
                if combo_id not in seen:
                    print combo_id
                seen[combo_id] = True



if __name__ == '__main__':
	main()

