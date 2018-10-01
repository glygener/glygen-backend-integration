import os,sys
import json
import glob
from optparse import OptionParser
import commands

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-o","--organismType",action="store",dest="organismType",help="The type of organism (human, mouse)")
    (options,args) = parser.parse_args()
    for file in ([options.organismType]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    organism = options.organismType
    config_json = json.loads(open("../conf/config-1.json", "r").read())
    in_file = config_json["pathinfo"]["glytoucan_input"]+organism+"_gtc6.txt"
    in_dir = config_json["pathinfo"]["glytoucan_images_input"]
    out_dir = config_json["pathinfo"]["glytoucan_images_output"]+organism+"_glycan_images"

    temp =  os.path.exists(out_dir)
    if temp == False:
        cmd = "mkdir %s"  % (out_dir)
        x = commands.getoutput(cmd)

    row_count = 0
    with open(in_file) as f:
        for row in f:
            if row_count==0:
                row_count = 1
                continue
            else:
                cmd = "cp -n %s %s"  % (in_dir+row.strip()+".png",out_dir)
                x = commands.getoutput(cmd)

if __name__ == '__main__':
    main()
