import requests
import glob
import json

#url = "https://dsapi.glygen.org/dataset/upload"         ## production server
#url = "https://beta-dsapi.glygen.org/dataset/upload"    ## beta server
url = "https://dsapi.tst.glygen.org/dataset/upload"     ## tst server
#url = "https://dsapi.dev.glygen.org/dataset/upload"     ## dev server

file_list = glob.glob("tmp/glyco_sites_unicarbkb.csv")
for in_file in file_list:
    files = {"file": open(in_file, 'rb')}
    req_obj = {"format":"csv", "qctype":"glyco_site_unicarbkb", "dataversion":"1.12.1"}
    r = requests.post(url, files=files, data=req_obj, verify=False)
    out_file = ".".join(in_file.split(".")[:-1]) + ".json"
    with open(out_file, "w") as FW:
        FW.write("%s\n" % (json.dumps(json.loads(r.text), indent=4)))
    print ("Saved response in %s" % (out_file))


