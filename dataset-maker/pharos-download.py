import os
import csv
import sys
import json
import commands
import glob
import gzip
import requests


import libgly





def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
 
    global path_obj
    path_obj = config_obj["pathinfo"]


    top = 1000
    for i in xrange(1, 10000):
        skip = (i-1)*1000
        url = "https://pharos.nih.gov/idg/api/v1/targets?top=%s&skip=%s" % (top, skip)
        out_file = path_obj["downloads"] + "pharos/current/batch.%s.json" % (i)
        res = requests.get(url)
        
        if res.content.strip() != "":
            res_json = json.loads(res.content)
            if "count" in res_json:
                if res_json["count"] == 0:
                    print "batch %s count=0, terminating!" % (i)
                    sys.exit()
                else:
                    with open(out_file, 'w') as FW:
                        print "downloaded batch %s" % (i)
                        FW.write("%s\n" % (res.content))
            else:
                print "batch %s bad content, terminating!" % (i)
                sys.exit()
        else:
            print "batch %s returned empty content, terminating!" % (i)
            sys.exit()


    return


                


if __name__ == '__main__':
        main()
