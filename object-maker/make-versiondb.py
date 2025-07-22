import os,sys
import string
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
import datetime
import pytz
import csvutil


__version__="1.0"
__status__ = "Dev"



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-v","--ver",action="store",dest="ver",help="2.0.2")

    (options,args) = parser.parse_args()
    for file in ([options.ver]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    global wrk_dir
    
    ver = options.ver
    
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    jsondb_dir = wrk_dir + "/jsondb/"

    ts = datetime.datetime.now(pytz.timezone('US/Eastern'))
    ts = ts.strftime('%Y-%m-%d %H:%M:%S %Z%z')
    doc = { "component":"data", "version":ver, "release_date":ts}
                           

    out_file = jsondb_dir + "/versiondb/data.json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(doc, indent=4)))

    
    log_file = "logs/make-versiondb.log"
    msg = "make-versiondb: final created: %s version objects" % (1)
    csvutil.write_log_msg(log_file, msg, "w")


            

if __name__ == '__main__':
	main()

