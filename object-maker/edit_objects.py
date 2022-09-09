import os,sys
import string
import subprocess
from optparse import OptionParser
import glob
import json
import libgly




def main():


    tr_dict = {"ln2releases":"ln2data/releases","ln2wwwdata":"ln2data/releases/data/REL"}
    in_file = "generated/misc/bco_changes_1.csv"
    changes_dict = {}
    data_frame = {}
    #libgly.load_sheet(data_frame, in_file, ",")
    #f_list = data_frame["fields"]
    #for row in data_frame["data"]:
    #    val_dict = {}
    #    for f in f_list:
    #        val_dict[f] = row[f_list.index(f)]
    #    changes_dict[val_dict["bco_id"]] = {
    #        "name":val_dict["new_bco_name"],
    #        "usability":val_dict["new_usability_domain"]
    #    }


    json_db = "bcodb"
    cmd = "ls -d /data/shared/glygen/releases/data/v-*/jsondb/%s/" % (json_db)
    d_list = subprocess.getoutput(cmd).split("\n")

    for d in d_list:
        rel = d.split("/")[-4]
        if os.path.isdir("tmp/%s/jsondb/%s" % (rel, json_db)) == False:
            x = subprocess.getoutput("mkdir -p tmp/%s/jsondb/%s" % (rel,json_db))
        rel_dir = "/data/shared/glygen/releases/data/%s/" % (rel)
        file_list = glob.glob(rel_dir + "jsondb/%s/*.json" % (json_db))
        for in_file in file_list:
            bco_id = in_file.split("/")[-1].split(".")[0]
            doc_str = open(in_file,"r").read()
            #for k in tr_dict:
            #    v = tr_dict[k].replace("REL", rel)
            #    doc_str = doc_str.replace(k, v)
            doc = json.loads(doc_str)
            #if bco_id in changes_dict:
            #    doc["provenance_domain"]["name"] = changes_dict[bco_id]["name"]
            #    doc["usability_domain"] = [changes_dict[bco_id]["usability"]]
            doc["object_id"] = "https://biocomputeobject.org/%s/%s" % (bco_id, rel)
            doc["provenance_domain"]["version"] = rel

            out_file = "tmp/%s/jsondb/%s/%s.json" % (rel, json_db, bco_id)            
            with open(out_file, "w") as FW:
                FW.write("%s\n" % (json.dumps(doc, indent=4)))
            #print ("Edited ... ", out_file)

        cmd = "mv tmp/%s/jsondb/%s/*.json %s" % (rel, json_db, d)
        x = subprocess.getoutput(cmd)
        print ("Edited and moved ... ", json_db, rel)

    





if __name__ == '__main__':
	main()



