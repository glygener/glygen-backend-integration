import os
import csv
import sys
import json
import commands
import glob
import gzip


import libgly



def remove_empty_files():

    cmd = "ls -l %spubchem/compound/current/sdf/*.sdf.gz*" % (path_obj["downloads"])
    lines = commands.getoutput(cmd).split("\n")
    for line in lines:
        file_size = line.split()[4]
        out_file = line.split()[-1]
        if file_size == "0":
            cmd = "rm %s" % (out_file)
            x = commands.getoutput(cmd)




def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
 
    global path_obj
    path_obj = config_obj["pathinfo"]


    if os.path.isdir("downloads/pubchem/compound/current/sdf/") == False:
        cmd = "mkdir -p downloads/pubchem/compound/current/sdf"
        x = commands.getoutput(cmd)
        cmd = "mkdir -p downloads/pubchem/compound/current/sdf4glygen"
        x = commands.getoutput(cmd)


    seen = {}
    url_tmplt = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
    url_tmplt += "Compound_%s_%s.sdf.gz"


    log_file = "logs/pubchem-downlads-step1.log"
    with open(log_file, "w") as FL:
        FL.write("started logging\n")

    #batch_size = 250000
    batch_size = 500000
    n_failed = 0
    for i in xrange(0, 1000):
        start = i*batch_size + 1
        end = start + batch_size - 1
        start = "000000000"[0:-len(str(start))] + str(start)
        end = "000000000"[0:-len(str(end))] + str(end)
        url = url_tmplt % (start, end)
        file_name = "Compound_%s_%s.sdf.gz" % (start, end)
        md5_file_name = "Compound_%s_%s.sdf.gz.md5" % (start, end)
        
        out_file_path = path_obj["downloads"] + "pubchem/compound/current/sdf/%s" % (file_name)
        md5_file_path = path_obj["downloads"] + "pubchem/compound/current/sdf/%s" % (md5_file_name)
        if os.path.isfile(out_file_path) == True:
            with open(log_file, "a") as FL:
                FL.write("%s exists\n" % (file_name))
        else:
            cmd = "wget -O %s %s" % (out_file_path, url)
            x = commands.getoutput(cmd)
            cmd = "wget -O %s %s" % (md5_file_path, url)
            x = commands.getoutput(cmd)
            with open(log_file, "a") as FL:
                FL.write("%s downloaded\n" % (file_name))
            remove_empty_files()



if __name__ == '__main__':
        main()
