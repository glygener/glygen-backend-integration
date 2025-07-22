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




def versionize_objects(ver):

    rel_dir = "/data/shared/glygen/releases/data/"
    
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

    versionize_objects(ver)

 
    wrk_dir = "/home/rykahsay/glygen-backend-integration/release-maker/"
    log_file = wrk_dir + "logs/frontend_release_progress.log"
   
    config_obj = json.loads(open("conf/config.json", "r").read()) 
    dir_dict = {
        "blastdb":"/data/projects/glygen/generated/datasets/blastdb/",
        "logs":"/data/projects/glygen/generated/datasets/logs/",
        "misc":"/data/projects/glygen/generated/misc/",
        "documentation":"/data/projects/glygen/documentation/",
        "glycanimages_snfg_png":"/data/projects/glygen/generated/datasets/reviewed/",
        "glycanimages_snfg_svg":"/data/projects/glygen/generated/datasets/reviewed/",
        "glycanimages_snfg_json":"/data/projects/glygen/generated/datasets/reviewed/"
    }
    
    for db in config_obj["frontend"]["dblist"]:
        k = "jsondb/%s" % (db)
        dir_dict[k] = "/data/projects/glygen/generated/datasets/jsondb/%s/"  % (db)


    file_list = glob.glob("/data/projects/glygen/generated/datasets/jsondb/motifdb/*.json")
    motif2gtc = {} 
    for in_file in file_list:
        doc = json.loads(open(in_file, "r").read())
        motif2gtc[doc["motif_ac"]] = doc["glytoucan_ac"]






    
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
        if subdir in ["glycanimages_snfg_png", "glycanimages_snfg_svg", "glycanimages_snfg_json"]:
            frmt = subdir.split("_")[-1]
            k = "glycanimages_snfg_" + frmt
            subdir_path = "%s/%s/" % (release_dir, subdir)

            src_tarball_file = "%s/glycan_images_snfg_extended_%s.tar.gz" % (dir_dict[subdir], frmt)
            cmd = "cp %s %s " % (src_tarball_file, subdir_path)
            x = subprocess.getoutput(cmd)
            #print (cmd)
            log_progress(log_file, "finished copying src: glycan_images_snfg_extended_%s.tar.gz" %(frmt), "a")
            
            os.chdir(release_dir + "/" + subdir)
            cmd = "pwd"
            x = subprocess.getoutput(cmd)
            #print (x)
            
            cmd = "rm -f *.%s" % (frmt)
            #print (cmd)
            x = subprocess.getoutput(cmd)
            log_progress(log_file, "finished removing glycan_images_snfg_extended_%s.tar.gz" % (frmt), "a")
 
            cmd = "tar xfz glycan_images_snfg_extended_%s.tar.gz" % (frmt)
            #print (cmd)
            x = subprocess.getoutput(cmd)
            log_progress(log_file, "finished untar: glycan_images_snfg_extended_%s.tar.gz" % (frmt) , "a")

            cmd = "rm -f glycan_images_snfg_extended_%s.tar.gz" % (frmt) 
            x = subprocess.getoutput(cmd)
            #print (cmd)
            log_progress(log_file, "finished removing glycan_images_snfg_extended_%s.tar.gz" % (frmt), "a")
            
            for motif_ac in motif2gtc:
                glytoucan_ac = motif2gtc[motif_ac]
                cmd = "ln -s %s.%s %s.%s" % (glytoucan_ac,frmt, motif_ac, frmt)
                x = subprocess.getoutput(cmd)

            log_progress(log_file, "finished ln for images_snfg_%s" % (frmt), "a") 

            os.chdir(wrk_dir)
        else:
            if os.path.isdir(subdir_path) == True:
                cmd = "rm -rf %s" % (subdir_path)
                x, y = subprocess.getstatusoutput(cmd)
            cmd = "cp -r %s %s" % (dir_dict[subdir], subdir_path)
            x, y = subprocess.getstatusoutput(cmd)
    
        nfiles = len(glob.glob("%s/*" % (subdir_path)))
        log_progress(log_file, "copied %s files to release dir: %s"% (nfiles, subdir), "a")


    cmd = "chmod -R 777 " + wrk_dir + "/logs"
    x, y = subprocess.getstatusoutput(cmd)



if __name__ == '__main__':
	main()

