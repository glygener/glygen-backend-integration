import os,sys
from optparse import OptionParser
import glob
import json
import subprocess
import datetime
import hashlib


__version__="1.0"
__status__ = "Dev"



def versionize_objects(ver):

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
        for obj in doc["description_domain"]["pipeline_steps"]:
            for o in obj["input_list"] + obj["output_list"]:
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
        

    return





def log_progress(log_file, log_msg, mode):

    with open(log_file, mode) as FW:
        ts = subprocess.getoutput("date")
        FW.write("... %s\n" % (log_msg))

    return

def get_file_dict(bcodb_dir):

    file_dict = {}
    doc_list = []
    for in_file in glob.glob(bcodb_dir + "/*.json"):
        doc_list.append(json.loads(open(in_file, "r").read()))

    for doc in doc_list:
        bco_id = doc["object_id"].split("/")[-2]
        tmp_list = []
        cn_doc = doc
        if "io_domain" in cn_doc:
            if "output_subdomain" in cn_doc["io_domain"]:
                if cn_doc["io_domain"]["output_subdomain"] != []:
                    for o in cn_doc["io_domain"]["output_subdomain"]:
                        file_name = o["uri"]["filename"].strip()
                        file_ext = file_name.split(".")[-1]
                        if file_name.find(".stat.csv") == -1 and file_ext != "log":
                            tmp_list.append(file_name)
        if tmp_list != []:
            file_dict[bco_id] = tmp_list

    return file_dict


###############################
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
    release_dir = "/data/shared/glygen/releases/data/v-%s/" % (ver)

    
    wrk_dir = "/home/rykahsay/glygen-backend-integration/release-maker/"
    log_file = wrk_dir + "logs/dsviewer_release_progress.log"
   
    config_obj = json.loads(open("conf/config.json", "r").read()) 

    dir_dict = {
        "reviewed":"/data/projects/glygen/generated/datasets/reviewed/", 
        "compiled":"/data/projects/glygen/generated/datasets/compiled/"
    }
    #dir_dict = {}
    #config_obj["dsviewer"]["dblist"] = ["recordsdb"] 
    for db in config_obj["dsviewer"]["dblist"]:
        k = "jsondb/%s" % (db)
        dir_dict[k] = "/data/projects/glygen/generated/datasets/jsondb/%s/"  % (db)

   
    bcodb_dir = "/data/projects/glygen/generated/datasets/jsondb/bcodb/"
    file_dict = get_file_dict(bcodb_dir)
    bco_id_list = list(file_dict.keys())
    
    log_progress(log_file, "creating directories ...", "w")
    if os.path.isdir(release_dir) == False:
        cmd = "mkdir  %s" % (release_dir)
        x, y = subprocess.getstatusoutput(cmd)
        log_progress(log_file, "finished creating main release folder", "a")

    for dir_name in dir_dict:
        dir_path = release_dir + "/%s/"  % (dir_name)
        if os.path.isdir(dir_path) == False:
            cmd = "mkdir -p %s" % (dir_path)
            x, y = subprocess.getstatusoutput(cmd)
    log_progress(log_file, "finished making release subfolders", "a")


    for subdir in dir_dict:
        subdir_path = release_dir + "/%s"  % (subdir)
        if os.path.isdir(subdir_path) == True:
            cmd = "rm -rf %s" % (subdir_path)
            x, y = subprocess.getstatusoutput(cmd)
        cmd = "cp -r %s %s" % (dir_dict[subdir], subdir_path)
        x, y = subprocess.getstatusoutput(cmd)
        nfiles = len(glob.glob("%s/*" % (subdir_path)))
        log_progress(log_file, "copied %s files to release dir: %s"% (nfiles, subdir), "a")


    versionize_objects(ver)
    log_progress(log_file, "finished versionizing objects", "a")

    


    cmd = "chmod -R 777 " + wrk_dir + "/logs"
    x, y = subprocess.getstatusoutput(cmd)





if __name__ == '__main__':
	main()

