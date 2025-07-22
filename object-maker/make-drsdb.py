#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import subprocess
from optparse import OptionParser
import libgly
import csvutil
import datetime
import pytz
import hashlib

__version__="1.0"
__status__ = "Dev"





def load_bcoid2fname(fname2bcoid, bcoid2fname):

    file_list = glob.glob("jsondb/bcodb/*.json")
    for in_file in file_list:
        bco_doc = json.loads(open(in_file, "r").read())
        bco_id = bco_doc["object_id"].split("/")[-2]
        desc = bco_doc["provenance_domain"]["name"]
        if "io_domain" in bco_doc:
            for obj in bco_doc["io_domain"]["output_subdomain"]:
                if "uri" in obj:
                    if "filename" in obj["uri"]:
                        file_name = obj["uri"]["filename"]
                        if file_name.strip() != "" and file_name.find(".stat.csv") == -1:
                            fname2bcoid[file_name] = bco_id
                            if bco_id not in bcoid2fname:
                                bcoid2fname[bco_id] = {}
                            bcoid2fname[bco_id][file_name] = desc
    return 



def main():


    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-v","--ver",action="store",dest="ver",help="2.0.2")

    (options,args) = parser.parse_args()
    for file in ([options.ver]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    
    ver = options.ver


    global wrk_dir
    global field_dict
    global io_dict 
    global generated_dir
    global log_file
 

    DEBUG = False
    #DEBUG = True

    generated_dir = "/data/projects/glygen/generated/"
    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"
    reviewed_dir = wrk_dir + "/reviewed/"
    log_file = "logs/make-drsdb.log"
    
    msg = "make-drsdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")
    
    fname2bcoid, bcoid2fname = {}, {}
    load_bcoid2fname(fname2bcoid, bcoid2fname)


    api_base_url = "https://dsapi.tst.glygen.org/drs/objects/"
    file_base_url = "https://data.glygen.org/ln2data/releases/data/v-%s/reviewed/" % (ver)
    ts = datetime.datetime.now(pytz.timezone('US/Eastern'))
    ts = ts.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    bco_id_list = sorted(list(bcoid2fname.keys()))
    #bco_id_list = bco_id_list[:2]
    


    for bco_id in bco_id_list:
        file_name_list = sorted(list(bcoid2fname[bco_id].keys()))
        for idx in range(0, len(file_name_list)):
            file_name = file_name_list[idx]
            file_desc = bcoid2fname[bco_id][file_name]
            file_ext = file_name.split(".")[-1]
            file_path = "reviewed/%s" % (file_name)
            if os.path.isfile(file_path) == False:
                continue
            drs_id = "%s_%s" % (bco_id, idx + 1)   
            mem_type = "text/csv"
            mem_type = "text/tsv" if file_ext in ["tsv"] else mem_type
            mem_type = "text/plain" if file_ext in ["fasta"] else mem_type
            file_size = os.path.getsize(file_path)
            file_url = "%s%s" % (file_base_url, file_name)
            self_uri = "%s%s" % (api_base_url, drs_id)
            md5 = hashlib.md5(open(file_path,'rb').read()).hexdigest()
            doc = {
                "id":drs_id, 
                "name":file_name, 
                "description": file_desc,
                "self_uri": self_uri,
                "size":file_size,
                "version":ver,
                "mime_type": mem_type,
                "created_time": ts, 
                "updated_time": ts,
                "checksums": [{"checksum": md5, "type": "md5"}],
                "access_methods": [
                    {
                        "type": "https", 
                        "access_url": {"url": file_url},
                        "access_id": drs_id
                    }
                ]
            }
            out_file = wrk_dir + "/jsondb/drsdb/%s.json" % (drs_id)
            with open(out_file,"w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            msg = "make-drsdb: created file: jsondb/drsdb/%s.json" % (drs_id)
            csvutil.write_log_msg(log_file, msg, "a")



    cmd = "chmod -R 775 " + wrk_dir + "/jsondb/drsdb/"
    x, y = subprocess.getstatusoutput(cmd)



if __name__ == '__main__':
    main()



