import os,sys
import string
from optparse import OptionParser
import glob
import json
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
    ftp_dir = "/data/shared/glyds/releases/ftp/v-%s/" % (ver)


    if os.path.isdir(ftp_dir) == True:
        cmd = "rm -rf %s" % (ftp_dir)
        x = subprocess.getoutput(cmd)
    
    cmd = "mkdir -p " + ftp_dir
    x = subprocess.getoutput(cmd)

    cmd = "cp " + wrk_dir + "generated/misc/readme_ftp.txt " + ftp_dir + "/README.txt"
    x = subprocess.getoutput(cmd)

    os.chdir(release_dir + "/jsondb/")
    d_list = glob.glob("*")
    os.chdir(release_dir)
    for d in d_list:
        if d[-2:] != "db":
            continue
        cmd = "tar cvf %s/%s.tar jsondb/%s" % (ftp_dir, d, d)
        x = subprocess.getoutput(cmd)
        cmd = "gzip %s/%s.tar" % (ftp_dir, d)
        x = subprocess.getoutput(cmd)
        print ("finished %s" % (d)) 

    d_list = ["reviewed", "compiled", "misc"]
    for d in d_list:
        cmd = "tar cvf %s/%s.tar %s" % (ftp_dir, d, d)
        x = subprocess.getoutput(cmd)
        cmd = "gzip %s/%s.tar" % (ftp_dir, d)
        x = subprocess.getoutput(cmd)
        print ("finished %s" % (d))







    


if __name__ == '__main__':
	main()

