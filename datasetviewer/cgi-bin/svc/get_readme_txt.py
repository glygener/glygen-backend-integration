import os,sys
import string
import commands
import csv
import traceback

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import json
import util

import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"




###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-o","--objid",action="store",dest="objid",help="Object ID")
    parser.add_option("-v","--objver",action="store",dest="objver",help="Object Version")


    (options,args) = parser.parse_args()
    for key in ([options.objid]):
        if not (key):
            parser.print_help()
            sys.exit(0)


    global config_json
    global db_obj
    global client


    out_json = {}

    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        db_obj = config_json[config_json["server"]]["dbinfo"]

    except Exception, e:
        out_json = {"taskstatus":0, "errormsg":"Loading config failed!"}

    try:
        client = MongoClient('mongodb://localhost:27017',
            username=db_obj["mongodbuser"],
            password=db_obj["mongodbpassword"],
            authSource=db_obj["mongodbname"],
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        dbh = client[db_obj["mongodbname"]]
        options.objid = options.objid.replace("GLYDS", "DSBCO_")
        obj_id = "http://data.glygen.org/%s" % (options.objid)
        obj_ver = options.objver
        query_obj = {"bco_id":obj_id}
        
        bco_collection = "c_bco_" + obj_ver
        doc = dbh[bco_collection].find_one(query_obj)

        print "<pre>"
        
        if doc == None:
            print "\n\n\tReadme file does not exist for object %s!" % (obj_id)
        else:
            creator_list = []
            author_list = []
            for o in doc["provenance_domain"]["contributors"]:
                if "createdBy" in o["contribution"]:
                    creator_list.append("%s [%s]" % (o["name"], o["email"]))
                if "contributedBy" in o["contribution"]:
                    author_list.append("%s [%s]" % (o["name"], o["email"]))

            print "<pre style=\"white-space: pre-wrap;\">"
            print "<b>Name</b>\n   - %s" % (doc["provenance_domain"]["name"])
            print "\n<b>Version</b>\n   - %s" % (doc["provenance_domain"]["version"])
            print "\n<b>Created on</b>\n  - %s" % (doc["provenance_domain"]["created"])
            print "\n<b>Modified on</b>\n   - %s" % (doc["provenance_domain"]["modified"])
            print "\n<b>Created by</b>\n   - %s" % (", ".join(creator_list))
            print "\n<b>Authors</b>\n   - %s" % (", ".join(author_list))
            print "\n<b>Checksum</b>\n   - %s" % (doc["checksum"])
            if doc["usability_domain"]  != []:
                print "\n<b>Usability Domain</b>\n   - %s" % ("\n   - ".join(doc["usability_domain"]))
            
            if doc["description_domain"]["keywords"] != [] or doc["description_domain"]["pipeline_steps"] != []:
                print "\n<b>Description Domain</b>"
                if "keywords" in doc["description_domain"]:
                    if doc["description_domain"]["keywords"] != []:
                        print "   Keywords\n   - %s" % ("\n   - ".join(doc["description_domain"]["keywords"]))
                if "pipeline_steps" in doc["description_domain"]:
                    if doc["description_domain"]["pipeline_steps"] != []:
                        print "\n   Pipeline Steps"
                        for o in doc["description_domain"]["pipeline_steps"]:
                            print "   - Step-%s: %s" % (o["step_number"], o["description"])


            if doc["execution_domain"] != {}:
                print "\n<b>Execution Domain</b>"
                if "scripts" in doc["execution_domain"]:
                    print "   Scripts"
                    for o in doc["execution_domain"]["script"]:
                        print "   - %s\n" % (o["uri"]["uri"])
                if "software_prerequisites" in doc["execution_domain"]:
                    print "   Software Prerequisites"
                    for o in doc["execution_domain"]["software_prerequisites"]:
                        print "      Name: %s" % (o["name"])
                        print "      Version: %s" % (o["version"])
                        print "      URI: %s" % (o["uri"]["uri"])
                        print ""

            if doc["io_domain"] != {}:
                print "<b>I/O Domain</b>"
                if "input_subdomain" in doc["io_domain"]:
                    if doc["io_domain"]["input_subdomain"] != []:
                        print "   Input Subdomain"
                        for o in doc["io_domain"]["input_subdomain"]:
                            print "      Name: %s" % (o["uri"]["filename"])
                            #print "      Media Type: %s" % (o["mediatype"])
                            print "      URI: %s" % (o["uri"]["uri"])
                            print ""
                if "output_subdomain" in doc["io_domain"]:
                    if doc["io_domain"]["output_subdomain"] != []:
                        print "   Output Subdomain"
                        for o in doc["io_domain"]["output_subdomain"]:
                            print "      Name: %s" % (o["uri"]["filename"])
                            #print "      Media Type: %s" % (o["mediatype"])
                            print "      URI: %s" % (o["uri"]["uri"])
                            print ""

            if "emperical_error" in doc["error_domain"]:
                if "outputfields" in doc["error_domain"]["emperical_error"]:
                    print "<b>Fields</b>"
                    for o in doc["error_domain"]["emperical_error"]["outputfields"]:
                        print "      Name: %s" % (o["fieldname"])
                        print "      Description: %s" % (o["fielddesc"])
                        print "      Unique Values: %s" % (o["uniquevalues"])
                        print ""   

            #print "    %s" % (json.dumps(doc["error_domain"]))
            #print ": %s" % ()
            #print ": %s" % ()
            #print ": %s" % ()
            #print ": %s" % ()
            print "</pre>"
        sys.exit()
    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

