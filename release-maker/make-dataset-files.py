import os,sys
import string
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient
import subprocess
import datetime


__version__="1.0"
__status__ = "Dev"




def log_progress(log_file, log_msg, mode):

    with open(log_file, mode) as FW:
        ts, y = subprocess.getstatusoutput("date")
        FW.write("... %s\t%s\n" % (ts, log_msg))

    return

def get_file_dict(jsondb_dir):

    file_dict = {}
    doc_list = []
    for in_file in glob.glob(jsondb_dir + "/bcodb/*.json"):
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
    log_file = "/software/glygen/logs/progress.log"
    dir_dict = {
        "reviewed":"/data/projects/glygen/generated/datasets/reviewed/", 
        "compiled":"/data/projects/glygen/generated/datasets/compiled/",
        "source":"/data/projects/glygen/generated/datasets/source/",
        "blastdb":"/data/projects/glygen/generated/datasets/blastdb/", 
        "jsondb":"/data/projects/glygen/generated/datasets/jsondb/", 
        "logs":"/data/projects/glygen/generated/datasets/logs/",  
        "misc":"/data/projects/glygen/generated/misc/", 
        "documentation":"/data/projects/glygen/documentation/",
        "glycanimages_snfg":"/data/projects/glygen/generated/datasets/reviewed/"
    }
    

    file_dict = get_file_dict(dir_dict["jsondb"])
    bco_id_list = list(file_dict.keys())
    


    rel_list = []
    bco_his = []
    
    #Create release directory and subdirectories
    log_progress(log_file, "creating directories ...", "w")
    if os.path.isdir(release_dir) == True:
        cmd = "rm -rf %s" % (release_dir)
        x, y = subprocess.getstatusoutput(cmd)
        log_progress(log_file, "finished clearing old release dir", "a")
    cmd = "mkdir  %s" % (release_dir)
    x, y = subprocess.getstatusoutput(cmd)
    for dir_name in dir_dict:
        dir_path = release_dir + "/%s/"  % (dir_name)
        cmd = "mkdir  %s" % (dir_path)
        x, y = subprocess.getstatusoutput(cmd)
    log_progress(log_file, "finished making release subfolders", "a")

    for dir_name in dir_dict:
        if dir_name in ["reviewed", "glycanimages_snfg"]:
            continue
        dir_path = release_dir + "/%s/"  % (dir_name)
        cmd = "cp -r %s %s" % (dir_dict[dir_name], release_dir)
        x, y = subprocess.getstatusoutput(cmd)
        log_progress(log_file, "finished copying %s" % (dir_dict[dir_name]), "a")
    

    #create glycanimages
    cmd = "cp %s/glycan_images_snfg_extended_png.tar.gz %s/glycanimages_snfg/ "
    cmd = cmd % (dir_dict["glycanimages_snfg"], release_dir)
    x, y = subprocess.getstatusoutput(cmd)
    log_progress(log_file, "finished copying glycan_images_snfg_extended_png.tar.gz", "a")

    os.chdir(release_dir + "/glycanimages_snfg/")
    cmd = "tar xfz glycan_images_snfg_extended_png.tar.gz"
    x, y = subprocess.getstatusoutput(cmd)
    log_progress(log_file, "finished untar glycan_images_snfg_extended_png.tar.gz", "a")
    
    cmd = "rm -f glycan_images_snfg_extended_png.tar.gz"
    x, y = subprocess.getstatusoutput(cmd)
    log_progress(log_file, "finished removing glycan_images_snfg_extended_png.tar.gz", "a")
    os.chdir(".")

    bco_count = len(bco_id_list)
    log_progress(log_file,"Copying reviewed files", "a")
    nfiles = 0
    for i in range(0, bco_count):
        bco_id = bco_id_list[i]
        for file_name in file_dict[bco_id]:
            ds_file = release_dir + "/reviewed/" + file_name
            file_name_prefix = ".".join(file_name.split(".")[0:-1])
            cmd = "cp %s/%s* %s/reviewed/" % (dir_dict["reviewed"], file_name_prefix, release_dir)
            x, y = subprocess.getstatusoutput(cmd)
            nfiles += 1

    log_progress(log_file, "copied %s reviewed files " %(nfiles), "a") 

    #Make release note 
    log_msg = "Creating release note file ..."
    log_progress(log_file, log_msg, "a")
    release_date = datetime.datetime.now().strftime('%m/%d/%Y')
    release_file = "%s/reviewed/release-notes.txt" % (release_dir)
    with open(release_file, "w") as FW:
        FW.write("v-%s %s\n" % (ver, release_date))

    os.chdir(".")
    for db_dir in glob.glob("%s/*" % (release_dir)):
        d = db_dir.split("/")[-1]
        if d in ["jsondb"]:
            continue
        nfiles = len(glob.glob("%s/*" % (db_dir)))
        log_progress(log_file, "copied %s files to %s/"% (nfiles, db_dir), "a")

    for d_one in glob.glob("%s/jsondb*" % (release_dir)):
        for db_dir in glob.glob("%s/*" % (d_one)):
            nfiles = len(glob.glob("%s/*.json" % ( db_dir)))
            log_progress(log_file, "copied %s files to %s/"% (nfiles, db_dir), "a")

    cmd = "chmod -R 777 /software/glygen/logs/"
    x, y = subprocess.getstatusoutput(cmd)



if __name__ == '__main__':
	main()

