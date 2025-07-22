#!/usr/bin/python
import os,sys
import string
import csv
import json
import glob
import subprocess
from optparse import OptionParser
import csvutil


__version__="1.0"
__status__ = "Dev"




def get_readme(bco_id, doc):

    if doc == None:
        out_buffer += "\n\n\tReadme file does not exist for object %s!" % (bco_id)
    else:
        creator_list = []
        author_list = []
        for o in doc["provenance_domain"]["contributors"]:
            for p in ["createdBy"]:
                if "contribution" in o:
                    if p in o["contribution"]:
                        email = o["email"] if "email" in o else ""
                        name = o["name"] if "name" in o else ""
                        creator_list.append("%s [%s]" % (name, email))
            for p in ["contributedBy", "authoredBy"]:
                if "contribution" in o:
                    if p in o["contribution"]:
                        email = o["email"] if "email" in o else ""
                        name = o["name"] if "name" in o else ""
                        author_list.append("%s [%s]" % (name, email))

        modified = ""
        if "modified" in doc["provenance_domain"]:
            modified = doc["provenance_domain"]["modified"]
        checksum = doc["checksum"] if "checksum" in doc else ""

        out_buffer = ""
        out_buffer += "Name\n   - %s" % (doc["provenance_domain"]["name"])
        out_buffer += "\n\nVersion\n   - %s" % (doc["provenance_domain"]["version"])
        if "created" in doc["provenance_domain"]:
            out_buffer += "\n\nCreated on\n  - %s" % (doc["provenance_domain"]["created"])
        out_buffer += "\n\nModified on\n   - %s" % (modified)
        out_buffer += "\n\nCreated by\n   - %s" % (", ".join(creator_list))
        out_buffer += "\n\nAuthors\n   - %s" % (", ".join(author_list))
        out_buffer += "\n\nChecksum\n   - %s" % (checksum)
        if doc["usability_domain"]  != []:
            out_buffer += "\n\nUsability Domain\n   - %s" % ("\n   - ".join(doc["usability_domain"]))

        if "description_domain" in doc:
            if "keywords" not in doc["description_domain"]:
                doc["description_domain"]["keywords"] = []
            if doc["description_domain"]["keywords"] != [] or doc["description_domain"]["pipeline_steps"] != []:
                out_buffer += "\n\nDescription Domain"
                if "keywords" in doc["description_domain"]:
                    if doc["description_domain"]["keywords"] != []:
                        if None in doc["description_domain"]["keywords"]:
                            doc["description_domain"]["keywords"].remove(None)
                        out_buffer += "\n   Keywords\n   - %s" % ("\n   - ".join(doc["description_domain"]["keywords"]))
                if "pipeline_steps" in doc["description_domain"]:
                    if doc["description_domain"]["pipeline_steps"] != []:
                        out_buffer += "\n\n   Pipeline Steps"
                        for o in doc["description_domain"]["pipeline_steps"]:
                            if "step_number" in o and "description" in o:
                                out_buffer += "\n   - Step-%s: %s" % (o["step_number"], o["description"])

        if doc["execution_domain"] != {}:
            out_buffer += "\n\nExecution Domain"
            if "script" in doc["execution_domain"]:
                out_buffer += "\n   Scripts"
                for o in doc["execution_domain"]["script"]:
                    out_buffer += "   - %s\n" % (o["uri"]["uri"])
            if "software_prerequisites" in doc["execution_domain"]:
                out_buffer += "\n\n   Software Prerequisites"
                for o in doc["execution_domain"]["software_prerequisites"]:
                    out_buffer += "\n      Name: %s" % (o["name"])
                    out_buffer += "\n      Version: %s" % (o["version"])
                    out_buffer += "\n      URI: %s" % (o["uri"]["uri"])
                    out_buffer += ""

        if doc["io_domain"] != {}:
            out_buffer += "\n\nI/O Domain"
            if "input_subdomain" in doc["io_domain"]:
                if doc["io_domain"]["input_subdomain"] != []:
                    out_buffer += "\n   Input Subdomain"
                    for o in doc["io_domain"]["input_subdomain"]:
                        if "uri" not in o:
                            continue
                        if "filename" not in o["uri"]:
                            continue
                        out_buffer += "\n      Name: %s" % (o["uri"]["filename"])
                        out_buffer += "\n      URI: %s" % (o["uri"]["uri"])
                        out_buffer += ""
            if "output_subdomain" in doc["io_domain"]:
                if doc["io_domain"]["output_subdomain"] != []:
                    out_buffer += "\n   Output Subdomain"
                    for o in doc["io_domain"]["output_subdomain"]:
                        out_buffer += "\n      Name: %s" % (o["uri"]["filename"])
                        out_buffer += "\n      URI: %s" % (o["uri"]["uri"])
                        out_buffer += ""
                        if o["uri"]["filename"].find(".stat.") != -1:
                            stat_file = reviewed_dir + "%s"
                            stat_file = stat_file % (o["uri"]["filename"])
                            if os.path.isfile(stat_file):
                                stat_lines = open(stat_file, "r").read().split("\n")
                                out_buffer += "\n\nFields"
                                for stat_line in stat_lines[1:]:
                                    if stat_line.strip() == "":
                                        continue
                                    s_p = stat_line.split(",")
                                    out_buffer += "\n      Name: %s" % (s_p[1])
                                    out_buffer += "\n     Description: %s" % (s_p[2])

                                    out_buffer += "\n      Unique Values: %s" % (s_p[0])
                                    out_buffer += ""


    return out_buffer


def get_preview(doc):

    out_json = []
    file_name = doc["io_domain"]["output_subdomain"][0]["uri"]["filename"].strip()
    file_type = file_name.split(".")[-1]
    file_path = reviewed_dir +  file_name
    if os.path.exists(file_path) == False:
        return {"type":"", "data":""}
    else:
        species_short = file_name.split("_")[0]
        if file_type in ["csv", "tsv"]:
            delim = "," if file_type == "csv" else "\t"
            out_json = []
            with open(file_path, 'r') as FR:
                csvGrid = csv.reader(FR, delimiter=delim, quotechar='"')
                rowCount = 0
                for row in csvGrid:
                    rowCount += 1
                    if rowCount == 1:
                        tmp_list = []
                        for val in row:
                            tmp_list.append({"type":"string", "label":val})
                        out_json.append(tmp_list)
                    else:
                        tmp_list = []
                        for val in row:
                            tmp_list.append(val)
                        out_json.append(tmp_list)
                    if rowCount == 3:
                        break
            return {"type":"table", "data":out_json}
        elif file_type == "fasta":
            html_cn = "<pre>"
            seqCount = 0
            in_file = reviewed_dir + file_name
            seq_count = 0
            with open(in_file, "r") as FR:
                for line in FR:
                    if line[0] == ">":
                        html_cn += "\n"
                        seq_count += 1
                        if seq_count == 10:
                            break
                    html_cn += line
            html_cn += "</pre>"
            return {"type":"html", "data":html_cn}
        elif file_type.lower() in ["gif", "jpeg", "jpg"]:
            url = release_root + "reviewed/%s" % (file_name)
            html_cn = "<img src=\"%s\"><br>" % (url)
            return {"type":"html", "data":html_cn}
        elif file_type.lower() in ["mp4"]:
            url = release_root + "reviewed/%s" % (file_name)
            html_cn = "<video controls=\"controls\"><source src=\"%s\" type=\"video/mp4\"></video><br>" % (url)
            return {"type":"html", "data":html_cn}
        elif file_type in ["gz"] and file_name.find("glycan_images") != -1:
            #xxxxx
            file_list = glob.glob(glycan_images_dir + "*.png")
            html_cn = ""
            for f in file_list[1000:1010]:
                recordId = f.split("/")[-1].split(".")[0]
                url = release_root + "/glycanimages_snfg/" +  f.split("/")[-1]
                html_cn +=  "%s <br>" % (recordId)
                html_cn += "<img src=\"%s\"><hr><br>" % (url) 
            return {"type":"html", "data":html_cn}
        elif file_type in ["gp", "gb", "nt"]:
            in_file = reviewed_dir + file_name
            html_cn = ""
            line_count = 0
            with open(in_file, "r") as FR:
                for line in FR:
                    html_cn += line 
                    line_count += 1
                    if line_count == 100:
                        break
            return {"type":"text", "data":html_cn}
        else:
            html_cn = "<pre>Please implement service for %s preview!</pre>" % (file_type)
            return {"type":"html", "data":html_cn}


    return out_json



def process_bco_docs(doc_list):

    new_doc_list = []
    for doc_obj in doc_list:
        #doc = doc_obj["fields"]["contents"]
        doc = doc_obj
        category_dict = {}
        if "extension_domain" in doc:
            for oo in doc["extension_domain"]:
                if "dataset_extension" in oo:
                    if "dataset_categories" in oo["dataset_extension"]:
                        for o in oo["dataset_extension"]["dataset_categories"]:
                            cat_name = o["category_name"].replace(" ", "_")
                            cat_value = o["category_value"].strip()
                            if cat_value.strip() == "":
                                continue
                            if cat_name not in category_dict:
                                category_dict[cat_name] = []
                            if cat_value not in category_dict[cat_name]:
                                category_dict[cat_name].append(cat_value)

        retired_flag = ""
        for s in ["dataset_status", "status"]:
            if s in category_dict:
                for v in category_dict[s]:
                    if v.lower() in ["retired"]:
                        retired_flag = "retired"
        if retired_flag == "retired":
            continue

        file_name_list = []
        if "io_domain" in doc:
            if "output_subdomain" in doc["io_domain"]:
                if doc["io_domain"]["output_subdomain"] != []:
                    for ooo in doc["io_domain"]["output_subdomain"]:
                        if "uri" in ooo:
                            if "filename" in ooo["uri"]:
                                file_name = ooo["uri"]["filename"].strip()
                                if file_name != "":
                                    if file_name.find(".stat.csv") == -1:
                                        file_name_list.append(file_name)
        if len(file_name_list) == 0:
            continue
        bco_id = doc_obj["object_id"].split("/")[-2]
        file_type = file_name_list[0].split(".")[-1]
        if category_dict == {}:
            category_dict["file_type"] = file_type.lower()
        if "usability_domain" not in doc:
            doc["usability_domain"] = []
        desc = doc["usability_domain"][0] if doc["usability_domain"] != [] else ""
        max_desc_len = 100
        if len(desc) > max_desc_len:
            word_list = desc.split(" ")
            desc = ""
            for word in word_list:
                if len(desc) > max_desc_len:
                    desc += "..."
                    break
                desc += "%s " % (word)

        obj = {
            "bcoid": bco_id,
            "title":doc["provenance_domain"]["name"],
            "description":desc,
            "filename_list":", ".join(file_name_list),
            "filename": file_name_list[0],
            "filetype": file_type,
            "categories":category_dict
        }
        if file_type in ["csv", "tsv", "txt"]:
            obj["minitable"] = get_mini_table(file_name_list[0], file_type)
        else:
            obj["iconfilename"] = "%s_icon_all.png" % (file_type)
      
        doc_obj["extract"] = obj

    return  


def get_mini_table(file_name, file_type):

    header_row, body_rows = [], []
    delim = "," if file_type == "csv" else "\t"
    file_path = reviewed_dir + file_name

    if os.path.isfile(file_path):
        with open(file_path, 'r') as FR:
            data_frame = csv.reader(FR, delimiter=delim, quotechar='"')
            row_count = 0
            for row in data_frame:
                row_count += 1
                if row == []:
                    continue
                if row_count == 1:
                    if len(row) == 1:
                        header_row = [row[0][0:12] + "..."]
                    else:
                        header_row = [row[0][0:12] + "...", row[1][0:12] + "..."]
                elif row_count < 4:
                    if len(row) == 1:
                        body_rows.append([row[0][0:12]])
                    else:
                        body_rows.append([row[0][0:12] + "...", row[1][0:12] + "..."])
                else:
                    break
    return {"content":body_rows, "headers":header_row, "colwidth": ["40%", "20%"]}





def main():


    global release_root
    global data_ver
    global reviewed_dir
    global glycan_images_dir
    global wrk_dir

    wrk_dir = "/data/shared/repos/glygen-backend-integration/object-maker/"

    data_ver = "x.x.x"
    reviewed_dir = wrk_dir + "/reviewed/"
    release_root = "/ln2data/releases/data/v-%s/" % (data_ver)
    glycan_images_dir = wrk_dir + "/downloads/glytoucan/current/export/snfg/extended/png/"


         
    log_file = "logs/make-extractdb.log"
    msg = "make-extractdb: started logging"
    csvutil.write_log_msg(log_file, msg, "w")

    doc_list = []
    file_list = glob.glob(wrk_dir + "/jsondb/bcodb/*.json")
    for in_file in sorted(file_list):
        file_name = in_file.split("/")[-1]
        bco_id = file_name.split(".json")[0]
        doc_list.append(json.loads(open(in_file, "r").read()))

    process_bco_docs(doc_list)

    extractdb_dir = wrk_dir + "/jsondb/extractdb/" 
    if os.path.isdir(extractdb_dir) == False:
        cmd = "mkdir " + extractdb_dir
        x, y = subprocess.getstatusoutput(cmd)
    cmd = "rm -f %s/*" % (extractdb_dir)
    x, y = subprocess.getstatusoutput(cmd)

    for doc in doc_list:
        bco_id = doc["extract"]["bcoid"]
        out_file = extractdb_dir + "/%s.json" % (bco_id)
        doc["extract"]["sampledata"] = get_preview(doc)
        doc["extract"]["readme"] = get_readme(bco_id, doc)
        doc["extract"]["downloadurl"] = release_root + "reviewed/" +  doc["extract"]["filename"]
        with open(out_file,"w") as FW:
            FW.write("%s\n" % (json.dumps(doc["extract"], indent=4)))
        msg = "make-extractdb: created jsondb/extractdb/%s.json" % (bco_id)
        csvutil.write_log_msg(log_file, msg, "a")

    cmd = "chmod -R 775 " + extractdb_dir
    x, y = subprocess.getstatusoutput(cmd)


if __name__ == '__main__':
    main()


