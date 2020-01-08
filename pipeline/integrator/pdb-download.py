import os
from optparse import OptionParser
import sys,csv
import urllib,urllib2
import collections
import json
import commands
import glob

sys.path.append('../../glytools/')
import libgly




__version__="1.0"
__status__ = "Dev"

#################################
def error_handler(url):

    try:
        request = urllib2.Request(url)
        contact = "" # Please set your email address here to help us debug in case of problems.
        request.add_header('User-Agent', 'Python %s' % contact)
        response = urllib2.urlopen(request)
        page = response.read()
        return page
    except:
        pass
        return 1
###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-f","--format_type",action="store",dest="format_type",help="The format of file (pdb,cif,fasta)")
    (options,args) = parser.parse_args()
    for file in ([options.format_type]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]
 
    global path_obj
    path_obj = config_obj["pathinfo"]


    format_type = options.format_type+"_format"


    pdbid_list = []
    for f in glob.glob(path_obj["downloads"] + "pdb/current/*-linked-*"):
        if f.split("-")[-1] == "details":
            continue
        with open(f, 'r') as FR:
            for line in FR:
                pdb_id = line.strip().lower()
                if pdb_id not in pdbid_list:
                    pdbid_list.append(pdb_id)

    count=0
    out_dir = path_obj["downloads"] +  "pdb/current/%s/" % (format_type)
    if os.path.isdir(out_dir) == False:
        cmd = "mkdir %s" % (out_dir)
        x = commands.getoutput(cmd)

    for id in pdbid_list:
        if format_type=="pdb_format":
            url = "https://files.rcsb.org/download/"+id+".pdb"
            out_file = out_dir+id+".pdb" 
        elif format_type=="cif_format":
            url = "https://files.rcsb.org/download/"+id+".cif"
            out_file = out_dir+id+".cif"
        elif format_type=="fasta_format":
            url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList="+id+"&compressionType=uncompressed"
            out_file = out_dir+id+".fasta"
        else:
            print "This script does not support this format"
            break
        page = error_handler(url)
        if page == 1:
            print id
            count+=1
            #page = error_handler(url)    
        else:
            with open(out_file,"wb") as pdbfile:
                pdbfile.write(page)

    print count

if __name__ == '__main__':
        main() 
