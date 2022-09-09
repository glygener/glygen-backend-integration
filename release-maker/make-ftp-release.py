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


def get_dsfile2object_map_rows():

    row_list = []
    row = ["object_type", "path_in_obj", "dataset_file_pattern"]
    row_list.append(row)

    for in_file in glob.glob(wrk_dir + "/generated/misc/*_sectioninfo.json"):
        obj = json.loads(open(in_file, "r").read())
        obj_type = in_file.split("/")[-1].split("_")[0]

        for path in obj:
            p_list = []
            for p in obj[path]["sheetlist"]:
                mol = p.split("_")[0]
                pp = "%s*.csv" % (p) if mol == "glycan" else "*_%s*.csv" % (p)
                p_list.append(pp.replace("__", "_"))
            row = [obj_type, path, "; ".join(p_list)]
            row_list.append(row)


    return row_list


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


    global wrk_dir

    wrk_dir = "/home/rykahsay/glygen-backend-integration/release-maker/"

    release_dir = "/data/shared/glygen/releases/data/v-%s/" % (ver)
    #ftp_dir = "/data/shared/glygen/releases/ftp/v-%s/" % (ver)
    ftp_dir = "/data/shared/glygen/releases/ftp/v-%s/" % ("x.x")


    if os.path.isdir(ftp_dir) == True:
        cmd = "rm -rf %s" % (ftp_dir)
        x = subprocess.getoutput(cmd)

    cmd = "mkdir  -p " + ftp_dir + "/dataset_files/"
    x = subprocess.getoutput(cmd)

    cmd = "cp " + wrk_dir + "generated/misc/readme_ftp.txt " + ftp_dir + "/README.txt"
    x = subprocess.getoutput(cmd)

    
    row_list = get_dsfile2object_map_rows()
    out_file = ftp_dir + "/dataset2object_map.csv"
    with open(out_file, "w") as FW: 
        for row in row_list:
            FW.write("\"%s\"\n" % ("\",\"".join(row)))


    file_list = glob.glob(release_dir + "/reviewed/*.*")
    for path in file_list:
        ignore = False
        for k in ["GLY_0",".stat.csv"]:
            if path.find(k) != -1:
                ignore = True
        if path[-3:] == ".nt":
            ignore = True
        if ignore == True:
            continue
        cmd = "cp %s %s/dataset_files/" % (path, ftp_dir)
        x = subprocess.getoutput(cmd)
    

    os.chdir(ftp_dir)
    cmd = "tar cvf dataset_files.tar dataset_files"
    x = subprocess.getoutput(cmd)
    cmd = "gzip dataset_files.tar"
    x = subprocess.getoutput(cmd)
    cmd = "rm -rf dataset_files"
    x = subprocess.getoutput(cmd)


    d_list = [
        "bcodb", "glycandb", "motifdb", "speciesdb", "diseasedb", 
        "proteindb", "alignmentdb","networkdb","publicationdb","sitedb"
    ]
    os.chdir(release_dir + "/jsondb/")
    for d in d_list:
        cmd = "tar cvf %s/%s.tar %s" % (ftp_dir,d, d)
        x = subprocess.getoutput(cmd)
        cmd = "gzip %s/%s.tar" % (ftp_dir,d)
        x = subprocess.getoutput(cmd)



    


if __name__ == '__main__':
	main()

