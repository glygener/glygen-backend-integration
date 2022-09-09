import os,sys
import string
import subprocess
from optparse import OptionParser
import glob
import json
from Bio import SeqIO
import gzip

import csvutil
import datetime



__version__="1.0"
__status__ = "Dev"


def get_file_dict(jsondb_dir):

    file_dict = {}
    doc_list = []
    for in_file in glob.glob(jsondb_dir + "/bcodb/*.json"):
        doc_list.append(json.loads(open(in_file, "r").read()))

    for doc in doc_list:
        bco_id = doc["object_id"].split("/")[-2]
        tmp_list = []
        cn_doc = doc
        if "io_domain" in cn_doc:
            if "output_subdomain" in cn_doc["io_domain"]:
                if cn_doc["io_domain"]["output_subdomain"] != []:
                    for o in cn_doc["io_domain"]["output_subdomain"]:
                        file_name = o["uri"]["filename"].strip()
                        file_ext = file_name.split(".")[-1]
                        if file_name.find(".stat.csv") == -1 and file_name.find(".stats.csv") == -1 and file_ext != "log":
                            tmp_list.append(file_name)
        if tmp_list == []:
            continue

        retired_flag = "active"
        if "extension_domain" in cn_doc:
            if "dataset_categories" in cn_doc["extension_domain"]:
                for o in cn_doc["extension_domain"]["dataset_categories"]:
                    if o == None:
                        continue
                    cat_name = o["category_name"].replace(" ", "_")
                    cat_value = o["category_value"].strip()
                    if cat_name.lower() in ["status", "dataset_status"] and cat_value.lower() == "retired":
                        retired_flag = "retired"
        if retired_flag != "active":
            continue
        file_dict[bco_id] = tmp_list


    return file_dict



def dump_record_counts(jsondb_dir, count_file, reviewed_dir):

    file_dict = get_file_dict(jsondb_dir)
    FW = open(count_file, "w")
    row = ["field_count", "row_count", "id_count", "bco_id", "file_name"]
    FW.write("\"%s\"\n" % ("\",\"".join(row)))
    for bco_id in file_dict:
        for file_name in file_dict[bco_id]:
            in_file = "%s/%s" % (reviewed_dir, file_name)
            file_format = in_file.split(".")[-1]
            if os.path.isfile(in_file) == True:
                field_count, row_count, id_count = get_record_count(in_file)
                row = [str(field_count),str(row_count), str(id_count), bco_id,file_name]
                FW.write("\"%s\"\n" % ("\",\"".join(row)))
    FW.close()

    return


def sort_release_list(tmp_list, reversed_flag):

    factor_list = [10000, 1000, 1]
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

    return release_list



def add_delta_list(bco_id, file_name_one, file_name_two, rel, prev_rel, out_json):

    in_file_one = data_release_dir + "v-%s/reviewed/%s" % (prev_rel, file_name_one)
    in_file_two = current_reviewed_dir + "/%s" % (file_name_two)

    if os.path.isfile(in_file_one) == False or os.path.isfile(in_file_two) == False:
        return
    file_ext_one = file_name_one.split(".")[-1].lower()
    file_ext_two = file_name_two.split(".")[-1].lower()

    if file_ext_one not in ["csv", "tsv"] or file_ext_two not in ["csv", "tsv"]:
        return


    data_frame_one = {}
    sep_one = "," if file_ext_one == "csv" else "\t"
    csvutil.load_sheet_as_dict(data_frame_one, in_file_one, sep_one, "")
    f_set_one = set(data_frame_one["fields"])
    record_set_one  = set(list(data_frame_one["data"].keys()))

    data_frame_two = {}
    sep_two = "," if file_ext_two == "csv" else "\t"
    csvutil.load_sheet_as_dict(data_frame_two, in_file_two, sep_two, "")
    f_set_two = set(data_frame_two["fields"])
    record_set_two  = set(list(data_frame_two["data"].keys()))

    out_json[bco_id][rel]["fields_added"] = list(f_set_two - f_set_one)
    out_json[bco_id][rel]["fields_removed"] = list(f_set_one - f_set_two)
    out_json[bco_id][rel]["ids_added"] = list(record_set_two - record_set_one)
    out_json[bco_id][rel]["ids_removed"] = list(record_set_one - record_set_two)

    return






def get_record_count(in_file):


    file_ext = in_file.split(".")[-1].lower()
    field_count, row_count, id_count = 1, 1, 1
    if file_ext in ["csv", "tsv"]:
        sep = "\t" if file_ext == "tsv" else ","
        field_count, row_count, id_count = csvutil.get_sheet_stats(in_file, sep)
    elif file_ext in ["fasta"]:
        id_count = len(list(SeqIO.parse(in_file, "fasta")))
    #elif file_ext == "nt":
        #with open(in_file, "r") as FR:
        #    for line in FR:
        #        n += 1
    #elif file_ext == "gz":
        #id_count = 1
        #with gzip.open(in_file, 'rb') as FR:
        #    for line in FR:
        #        n += 1
    #else:
    #    id_count = 1

    return field_count, row_count, id_count







def dump_bco_history(pair):

    out_json = {}
    in_file = historydb_dir + "/record_count.csv"
    out_dir = historydb_dir
    if pair[-1] != current_rel:
        out_dir = data_release_dir + "v-%s/jsondb/historydb/" % (pair[-1])
        in_file = data_release_dir + "v-%s/jsondb/historydb/record_count.csv" % (pair[-1])
    cmd = "rm -f %s/*.json" % (out_dir)
    subprocess.getoutput(cmd)

    current_bco_list = []

    data_frame = {}
    csvutil.load_sheet(data_frame, in_file, [], ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        bco_id = row[f_list.index("bco_id")]
        current_bco_list.append(bco_id)

    for i in range(0, len(pair)):
        rel = pair[i]
        in_file, rel_date = "", datetime.datetime.now().strftime('%m/%d/%Y') 
        if rel == current_rel:
            in_file = historydb_dir + "/record_count.csv"
        else:
            in_file = data_release_dir + "v-%s/jsondb/historydb/record_count.csv" % (rel)
            notes_file = data_release_dir + "v-%s/reviewed/release-notes.txt" % (rel)
            rel_date = subprocess.getoutput("cat " + notes_file).strip().split(" ")[-1]
        if os.path.isfile(in_file) == False:
            continue
        data_frame = {}
        csvutil.load_sheet(data_frame, in_file, [], ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            bco_id = row[f_list.index("bco_id")]
            if bco_id not in current_bco_list:
                continue
            file_name = row[f_list.index("file_name")]
            field_count = int(row[f_list.index("field_count")])
            row_count = int(row[f_list.index("row_count")])
            id_count = int(row[f_list.index("id_count")])
            if bco_id not in out_json:
                out_json[bco_id] = {}
            if rel not in out_json[bco_id]:
                out_json[bco_id][rel] = {"file_name":file_name, "field_count":field_count,
                        "row_count":row_count,"id_count":id_count,"release_date":rel_date}
    
    for bco_id in sorted(out_json):
        for i in range(1, len(pair)):
            prev_rel = pair[i-1]
            rel = pair[i]
            field_count, prev_field_count = 0, 0
            row_count, prev_row_count = 0, 0
            id_count, prev_id_count = 0, 0
            if rel in out_json[bco_id]:
                field_count = out_json[bco_id][rel]["field_count"]
                row_count = out_json[bco_id][rel]["row_count"]
                id_count = out_json[bco_id][rel]["id_count"]
            if prev_rel in out_json[bco_id]:
                prev_field_count = out_json[bco_id][prev_rel]["field_count"]
                prev_row_count = out_json[bco_id][prev_rel]["row_count"]
                prev_id_count = out_json[bco_id][prev_rel]["id_count"]
            if rel in out_json[bco_id]:
                out_json[bco_id][rel]["fields_added"] = []
                out_json[bco_id][rel]["fields_removed"] = []
                out_json[bco_id][rel]["ids_added"] = []
                out_json[bco_id][rel]["ids_removed"] = []

            if prev_id_count != 0 and rel in out_json[bco_id] and prev_rel in out_json[bco_id]:
                file_name_one = out_json[bco_id][prev_rel]["file_name"]
                file_name_two = out_json[bco_id][rel]["file_name"]
                add_delta_list(bco_id,file_name_one, file_name_two, rel, prev_rel, out_json)
            prev_rel = rel
        out_file = "%s%s.pairs.json" % (out_dir,bco_id)
        with open(out_file, "w") as FW:
            tmp_json = {"bcoid":bco_id, "doctype":"pairs", "history":{}}
            for rel in sorted(out_json[bco_id], reverse=True):
                rel_key = rel.replace(".", "_")
                tmp_json["history"][rel_key] = out_json[bco_id][rel]
            FW.write("%s\n" % (json.dumps(tmp_json, indent=4, sort_keys=True)))

    return



###############################
def main():


    global data_release_dir
    global current_reviewed_dir
    global current_rel
    global historydb_dir

    current_rel = "x.x.x"
    wrk_dir = "/home/rykahsay/glygen-backend-integration/object-maker"
    data_release_dir = "/data/shared/glygen/releases/data/"
    current_reviewed_dir = wrk_dir + "/reviewed/"
    
    jsondb_dir = wrk_dir +  "/jsondb/"
    historydb_dir = wrk_dir +  "/jsondb/historydb/"
    current_count_file = historydb_dir + "/record_count.csv"


    url = "https://api.glygen.org//misc/verlist/"
    cmd = "curl -s -k %s" % (url)
    x,res = subprocess.getstatusoutput(cmd)
    tmp_list = json.loads(res)
    release_list = sort_release_list(tmp_list, False)

    prev_rel = release_list[-1]
    pair = [prev_rel, current_rel]
    #print ("doing pair ", pair)
    #exit()

    
    cmd = "rm -f " + historydb_dir + "/*.json"
    x, y = subprocess.getstatusoutput(cmd)
    
    dump_record_counts(jsondb_dir, current_count_file, current_reviewed_dir)
    dump_bco_history(pair)






if __name__ == '__main__':
	main()

