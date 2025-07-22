import os,sys
import string
import subprocess
from optparse import OptionParser
import glob
import json
import requests
import csvutil



__version__="1.0"
__status__ = "Dev"


def sort_release_list(tmp_list, reversed_flag):

    factor_list = [100000000, 1000, 1]
    rel_dict = {}
    for rel in tmp_list:
        parts = rel.split(".")
        ordr = 0
        for i in range(0,len(parts)):
            ordr += factor_list[i]*int(parts[i])
        rel_dict[ordr] = rel
    
    release_list = []

    for ordr in sorted(rel_dict, reverse=reversed_flag):
        release_list.append(rel_dict[ordr])


    seen = {}
    filtered_release_list = []
    for r in release_list:
        rr = r.split(".")[0] + "." + r.split(".")[1]
        if rr not in seen:
            seen[rr] = True
            filtered_release_list.append(r)


    return filtered_release_list




def main():

    global jsondb_dir
    global wrk_dir 

    ver = "x.x.x"

    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker"
    jsondb_dir = wrk_dir + "/jsondb/"

    in_file = wrk_dir + "/conf/init.json"
    out_json = json.loads(open(in_file, "r").read())
    out_json["dataversion"] = ver
    
    
    #url = "https://api.glygen.org/misc/verlist"
    #res = requests.get(url,  verify=False)
    #res_obj = json.loads(res.content)
    #release_list = sort_release_list(res_obj, True)

    tmp_list = json.loads(open("conf/prd_versions.json", "r").read())
    release_list = sort_release_list(tmp_list, True)


    if ver not in release_list:
        release_list = [ver] + release_list
    idx = release_list.index(ver)
    release_list = release_list[idx:]


    out_json["versionlist"] = release_list

    if os.path.isdir(jsondb_dir + "/initdb") == False:
        cmd = "mkdir " + jsondb_dir + "/initdb"
        x, y = subprocess.getstatusoutput(cmd)


    out_file = jsondb_dir + "/initdb/init.json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(out_json, indent=4)))

    cmd = "chmod -R 775 " + jsondb_dir + "/jsondb/initdb"
    x, y = subprocess.getstatusoutput(cmd)

    
    log_file = "logs/make-initdb.log"
    msg = "make-initdb: ... final created: %s path objects" % (1)
    csvutil.write_log_msg(log_file, msg, "w")






if __name__ == '__main__':
	main()



