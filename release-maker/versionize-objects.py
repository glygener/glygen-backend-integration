import os,sys
import string
import subprocess
from optparse import OptionParser
import glob
import json
import requests
import hashlib





def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog version___")
    parser.add_option("-v","--ver",action="store",dest="ver",help="Version")

    (options,args) = parser.parse_args()
    for key in ([options.ver]):
        if not (key):
            parser.print_help()
            sys.exit(0)

    ver = options.ver

    spec_version = "1.3.0"

    rel_dir = "/data/shared/glygen/releases/data/"
    
    # versionize bcodb objects
    file_list = glob.glob(rel_dir + "v-%s/jsondb/bcodb/*.json" % (ver))
    for obj_file in file_list:
        doc = json.loads(open(obj_file, "r").read())
        doc["object_id"] = doc["object_id"].replace("/v-x.x.x", "/v-%s"%(ver))
        doc["spec_version"] = spec_version
        doc["provenance_domain"]["version"] = "v-%s" % (ver)
        for obj in doc["io_domain"]["output_subdomain"]:
            obj["uri"]["uri"] = obj["uri"]["uri"].replace("/v-x.x.x/", "/v-%s/"%(ver))
        for obj in doc["io_domain"]["input_subdomain"]:
            obj["uri"]["uri"] = obj["uri"]["uri"].replace("/v-x.x.x/", "/v-%s/"%(ver))
        if "description_domain" in doc:
            if "pipeline_steps" in doc["description_domain"]:
                for obj in doc["description_domain"]["pipeline_steps"]:
                    for o in obj["input_list"] + obj["output_list"]:
                        if o == None:
                            continue
                        if "uri" in o:
                            o["uri"] = o["uri"].replace("/v-x.x.x/", "/v-%s/"%(ver))
        
        doc["etag"] = ""
        hash_obj = hashlib.md5(json.dumps(doc).encode('utf-8'))
        doc["etag"] = hash_obj.hexdigest()
                        

        with open(obj_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))
   
    # versionize extractdb objects
    file_list = glob.glob(rel_dir + "v-%s/jsondb/extractdb/*.json" % (ver))
    for obj_file in file_list:
        doc = json.loads(open(obj_file, "r").read())
        for k in ["downloadurl", "readme"]:
            doc[k] = doc[k].replace("v-x.x.x", "v-%s"%(ver))
        with open(obj_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))


    # versionize hitorydb objects
    file_list = glob.glob(rel_dir + "v-%s/jsondb/historydb/*.pairs.json" % (ver))
    file_list += glob.glob(rel_dir + "v-%s/jsondb/historydb/*.track.json" % (ver))
    for obj_file in file_list:
        doc = json.loads(open(obj_file, "r").read())
        new_ver = ver.replace(".", "_")
        if "x_x_x" in doc["history"]:
            doc["history"][new_ver] = doc["history"]["x_x_x"]
            doc["history"].pop("x_x_x")
        with open(obj_file, "w") as FW:
            FW.write("%s\n" % (json.dumps(doc, indent=4)))


    # versionize initdb objects
    obj_file = rel_dir + "v-%s/jsondb/initdb/init.json" % (ver)
    doc = json.loads(open(obj_file, "r").read())
    doc["versionlist"][0] = ver
    doc["dataversion"] = ver
    with open(obj_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(doc, indent=4)))
        




if __name__ == '__main__':
	main()



