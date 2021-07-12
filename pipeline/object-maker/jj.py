import glob
import json

json_file = "junk.json"
doc = json.loads(open(json_file,"r").read())
for obj in doc["phosphorylation"]:
    print obj
exit()


seen = {}
file_list = glob.glob("jsondb/alignmentdb/*.json")
for json_file in file_list:
    doc = json.loads(open(json_file,"r").read())
    for obj in doc["sequences"]:
        if len(obj["uniprot_id"]) not in seen:
            print len(obj["uniprot_id"]), obj["uniprot_id"]
        seen[len(obj["uniprot_id"])] = True
