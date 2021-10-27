#!/usr/bin/python
import os,sys
import string
import commands
import csv
import traceback

from Bio import SeqIO
from Bio.Seq import Seq


from optparse import OptionParser
import glob
from bson import json_util, ObjectId
import cgi
import json
import util
import datetime
import pytz
from pytz import timezone

import pymongo
from pymongo import MongoClient



__version__="1.0"
__status__ = "Dev"

def get_field_list(tmp_file, ext):

    if ext.lower() in ["csv", "tsv"]:
        delim = "," if ext == "csv" else "\t"
        row_list = []
        with open(tmp_file, "r") as FR:
            data_frame = csv.reader(FR, delimiter=delim, quotechar='"')
            row_count = 0
            field_count = 0
            for row in data_frame:
                row_count += 1
                if row_count == 1:
                    return row
    
    return []



def perform_sitemapping_qc(tmp_file, ext, in_json, mapping_rows, preview_rows):

    seen = {"ac":{}, "mainid":{},  "isoform":{}, "userfield":{}}
    pattern = data_path + current_release + "/reviewed/*_protein_canonicalsequences.fasta"


    aa_dict = {}
    in_file = data_path + current_release + "/misc/aadict.csv"
    lines = open(in_file, "r").read().split("\n")[1:]
    for line in lines:
        if line.strip() == "":
            continue
        row = line.split(",")
        row[1] = row[1].replace("\"", "")
        row[2] = row[2].replace("\"", "")
        aa_dict[row[1]] = row[2]
        aa_dict[row[1].lower()] = row[2].upper()
        aa_dict[row[1].upper()] = row[2].upper()
        aa_dict[row[2].upper()] = row[2].upper()


    seq_dict = {}
    file_list = glob.glob(pattern)
    for in_file in file_list:
        for record in SeqIO.parse(in_file, "fasta"):
            canon = record.id.split("|")[1]
            ac = canon.split("-")[0]
            seq_dict[canon] = str(record.seq.upper())
            seq_dict[ac] = str(record.seq.upper())

    seen = {}
    #delim = "," if ext == "csv" else "\t"
    delim = ","

    with open(tmp_file, "r") as FR:
        data_frame = csv.reader(FR, delimiter=delim, quotechar='"')
        row_count = 0
        field_count = 0
        fail_count = 0
        f_list = []
        for row in data_frame:
            row_count += 1
            if row_count == 1:
                f_list = row
                preview_rows.append(row + ["qc_comment"])
                preview_rows.append(["string" for v in row] + ["string"])
            else:
                ac = row[f_list.index(in_json["idfield"])]
                pos = row[f_list.index(in_json["posfield"])]
                residue = row[f_list.index(in_json["residuefield"])]
                if ac not in seq_dict:
                    comment = "unmapped protein ID"
                    combo_id = "%s %s %s %s" % (ac, pos, residue, comment)
                    if combo_id not in seen:
                        mapping_rows.append([ac, pos, residue, comment])
                        seen[combo_id] = True
                    continue
                
                if residue.strip() == "":
                    comment = "empty residue"
                    combo_id = "%s %s %s %s" % (ac, pos, residue, comment)
                    if combo_id not in seen:
                        mapping_rows.append([ac, pos, residue, comment])
                        seen[combo_id] = True
                    preview_rows.append(row + [comment])
                    continue

                if residue.strip().upper() not in aa_dict:
                    comment = "non-standard residue"
                    combo_id = "%s %s %s %s" % (ac, pos, residue, comment)
                    if combo_id not in seen:
                        mapping_rows.append([ac, pos, residue, comment])
                        seen[combo_id] = True
                    preview_rows.append(row + [comment])
                    continue

                if pos.isdigit() == False:
                    comment = "bad protein position"
                    combo_id = "%s %s %s %s" % (ac, pos, residue, comment)
                    if combo_id not in seen:
                        mapping_rows.append([ac, pos, residue, comment])
                        seen[combo_id] = True
                    preview_rows.append(row + [comment])
                    continue

                if int(pos) < 1 or int(pos) > len(seq_dict[ac]):
                    comment = "position out of sequence range"
                    combo_id = "%s %s %s %s" % (ac, pos, residue, comment)
                    if combo_id not in seen:
                        mapping_rows.append([ac, pos, residue, comment])
                        seen[combo_id] = True
                    preview_rows.append(row + [comment])
                    continue

                if seq_dict[ac][int(pos)-1] != aa_dict[residue]:
                    comment = "residue mismatch"
                    combo_id = "%s %s %s %s" % (ac, pos, residue, comment)
                    if combo_id not in seen:
                        mapping_rows.append([ac, pos, residue, comment])
                        seen[combo_id] = True
                    preview_rows.append(row + [comment])
                    continue

                preview_rows.append(row + [""])
    return



def perform_idmapping_qc(tmp_file, ext, in_json, mapping_rows, preview_rows):

    seen = {"ac":{}, "mainid":{},  "isoform":{}, "userfield":{}}
    pattern = data_path + current_release + "/reviewed/*_protein_masterlist.csv"
    if in_json["fieldtype"] == "glycan":
        pattern = data_path + current_release + "/reviewed/glycan_masterlist.csv"

    main_field = "uniprotkb_canonical_ac"
    if in_json["fieldtype"] == "glycan":
        main_field = "glytoucan_ac"

    file_list = glob.glob(pattern)
    for in_file in file_list:
        with open(in_file, "r") as FR:
            data_frame = csv.reader(FR, delimiter=",", quotechar='"')
            row_count = 0
            for row in data_frame:
                row_count += 1
                if row_count == 1:
                    f_list = row
                else:
                    main_id = row[f_list.index(main_field)]
                    seen["mainid"][main_id] = True
                    if in_json["fieldtype"] == "protein":
                        ac = main_id.split("-")[0]
                        seen["ac"][ac] = True
                        for f in ["reviewed_isoforms", "unreviewed_isoforms"]:
                            isoform = row[f_list.index(f)]
                            if isoform != "":
                                seen["isoform"][isoform] = main_id

    
    #delim = "," if ext == "csv" else "\t"
    delim = ","

    with open(tmp_file, "r") as FR:
        data_frame = csv.reader(FR, delimiter=delim, quotechar='"')
        row_count = 0
        field_count = 0
        fail_count = 0
        f_list = []
        for row in data_frame:
            row_count += 1
            if row_count == 1:
                f_list = row
                preview_rows.append(row + ["qc_comment"])
                preview_rows.append(["string" for v in row] + ["string"])
            else:
                user_field = row[f_list.index(in_json["fieldname"])]
                if in_json["fieldtype"] == "protein":
                    comment = ""
                    if user_field not in seen["mainid"] and user_field not in seen["ac"]:
                        comment = "Not in GlyGen ID space"
                        if user_field in seen["isoform"]:
                            comment = ", is an isoform of %s" % (seen["isoform"][user_field])
                        if user_field not in seen["userfield"]:
                            mapping_rows.append([user_field, comment])
                            seen["userfield"][user_field] = True
                    preview_rows.append(row + [comment])
                elif in_json["fieldtype"] == "glycan":
                    comment = ""
                    if user_field not in seen["mainid"]:
                        comment = "Not in GlyGen ID space"
                        if user_field not in seen["userfield"]:
                            mapping_rows.append([user_field, ""])
                            seen["userfield"][user_field] = True
                    preview_rows.append(row + [comment])



    return



def perform_format_sanity_qc(tmp_file, ext, preview_rows, sanity_rows,qcreport_rows):
    delim = "," if ext == "csv" else "\t"

    file_size = str(os.path.getsize(tmp_file)) + " Bytes"
    if ext.lower() in ["csv", "tsv"]:
        with open(tmp_file, "r") as FR:
            data_frame = csv.reader(FR, delimiter=delim, quotechar='"')
            row_count = 0
            field_count = 0
            fail_count = 0
            for row in data_frame:
                row_count += 1
                if row_count == 1:
                    field_count = len(row)
                    preview_rows.append(row)
                    preview_rows.append(["string" for v in row])
                elif len(row) != field_count:
                    record_str =   ", ".join(row)
                    sanity_rows.append(["format_sanity", record_str])
                    fail_count += 1
                else:
                    preview_rows.append(row)
            qcreport_rows.append(["file_size", file_size])
            qcreport_rows.append(["rows", str(row_count)])
            qcreport_rows.append(["columns", str(field_count)])
            qcreport_rows.append(["failed_rows", str(fail_count)])
    elif ext.lower() in ["png", "jpeg", "jpg"]:
        qcreport_rows.append(["file_size", file_size])


    return 




###############################
def main():
    
    
    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    msg = "Input JSON text"
    parser.add_option("-j","--injson",action="store",dest="injson",help=msg)


    form_dict = cgi.FieldStorage()
    (options,args) = parser.parse_args()
    local_flag = False
    in_json = {}
    if len(form_dict.keys()) > 0:
        in_json = json.loads(form_dict["injson"].value) if "injson" in form_dict else {}
    else:
        local_flag = True
        for key in ([options.injson]):
            if not (key):
                parser.print_help()
                sys.exit(0)
        in_json = json.loads(options.injson)



    global config_json
    global db_obj
    global client
    global current_release
    global data_path
    global data_root

    print "Content-Type: text/html"
    print
    
    #out_buffer += json.dumps(in_json, indent=4)
    #sys.exit()


    out_json = {}

    try:
        config_json = json.loads(open("conf/config.json", "r").read())
        custom_config_json = json.loads(open("conf/config.custom.json", "r").read())
        db_obj = custom_config_json[config_json["server"]]["dbinfo"]
        path_obj = custom_config_json[config_json["server"]]["pathinfo"]
        root_obj = custom_config_json[config_json["server"]]["rootinfo"]
        current_release = "v-" + config_json["datarelease"]
        data_path = path_obj["htmlpath"] + "/ln2releases/"
        data_root = root_obj["htmlroot"] + "/ln2releases/"

    except Exception, e:
        out_json = {"taskstatus":0, "errormsg":"Loading config failed!"}


    try:
        DEBUG = False
        #DEBUG = True
        tmp_file, ext = "", ""
        if DEBUG:
            file_name = "example.csv"
            ext = file_name.split(".")[-1]
            tmp_file = path_obj["uploads"] + "tmp/" + file_name

        sanity_rows = []
        preview_rows = []
        qcreport_rows = []
        mapping_rows = []
        if "action" not in in_json:
            in_json["action"] = "perform_format_sanity"
            if DEBUG == False:
                file_obj = form_dict["userfile"]
                ts = datetime.datetime.now(pytz.timezone('US/Eastern'))
                ts = ts.strftime('%m-%d-%Y')
                ext = file_obj.filename.split(".")[-1]
                fname = ".".join(file_obj.filename.split(".")[:-1])
                file_name = "%s.%s.%s" % (fname, ts, ext)
                tmp_file = path_obj["uploads"] + "tmp/" + file_name
                with open(tmp_file, "w") as FW:
                    FW.write(file_obj.file.read())
            sanity_rows.append(["qc_name","data_record"])
            sanity_rows.append(["string","string"])
            qcreport_rows.append(["field_name","field_value"])
            qcreport_rows.append(["string","string"])
            
            perform_format_sanity_qc(tmp_file, ext, preview_rows, sanity_rows,qcreport_rows)
            field_list = get_field_list(tmp_file, ext)
            file_info = {"filename":file_name, "fieldlist":field_list}
            out_json = {"previewrows":preview_rows, "sanityrows":sanity_rows, 
                    "qcreportrows":qcreport_rows,    
                "fileinfo":file_info}
            if ext.lower() in ["csv", "tsv"]:
                #Over write tmp_file with additiona qc_comment
                with open(tmp_file, "w") as FW:
                    for row in [preview_rows[0]] + preview_rows[2:]:
                        FW.write("\"%s\"\n" % ("\",\"".join(row)))

        elif in_json["action"] == "perform_idmapping":
            mapping_rows.append(["unmapped_id","comment"])
            mapping_rows.append(["string","string"])
            tmp_file = path_obj["uploads"] + "tmp/" + in_json["filename"]
            ext = in_json["filename"].split(".")[-1]
            perform_idmapping_qc(tmp_file, ext, in_json,mapping_rows,preview_rows)
            field_list = get_field_list(tmp_file, ext)
            file_info = {"filename":in_json["filename"], "fieldlist":field_list}
            out_json = {"mappingrows":mapping_rows, "previewrows":preview_rows, 
                    "fileinfo":file_info}
            #Over write tmp_file with additiona qc_comment
            with open(tmp_file, "w") as FW:
                for row in [preview_rows[0]] + preview_rows[2:]: 
                    FW.write("\"%s\"\n" % ("\",\"".join(row)))
        elif in_json["action"] == "say_hello":
            app_dir = "/software/test/"
            cmd = "matlab -nodisplay -nosplash -nodesktop -r "
            cmd += " \"cd %s;hello;exit;\" " % (app_dir)
            #cmd = "python /software/test/test.py"
            out_json = {"app_dir":app_dir,"last_line":commands.getoutput(cmd).split("\n")[-1]}

        elif in_json["action"] == "glycan_finder":
            tmp_file = path_obj["uploads"] + "tmp/" + in_json["filename"]
            ext = in_json["filename"].split(".")[-1]
            file_info = {"filename":in_json["filename"]}
            app_dir = "/software/glycan_finder/"
            cmd = "matlab -nodisplay -nosplash -nodesktop -r "
            cmd += " \"cd %s;inFile='%s';ClassOne;exit;\" " % (app_dir,tmp_file)
            glycan_list = commands.getoutput(cmd).split("\n")[-1].split(",")
            out_json = {
                "fileinfo":file_info, 
                "mappingrows":[
                    ["GlyToucan Accession", "Glycan Image"],
                    ["string", "string"]
                ]
            }
            for ac in glycan_list:
                link_one = "<a href=\"https://glygen.org/glycan/%s\" target=_>%s</a>" % (ac,ac)
                link_two = "<a href=\"https://gnome.glyomics.org/restrictions/GlyGen.StructureBrowser.html?focus=%s\" target=_>related glycans</a>" % (ac)
                img = "<img src=\"https://api.glygen.org/glycan/image/%s\">" % (ac)
                links = "%s (other %s)" % (link_one, link_two)
                out_json["mappingrows"].append([links, img])

        elif in_json["action"] == "perform_sitemapping":
            mapping_rows.append(["protein_id", "position", "residue", "comment"])
            mapping_rows.append(["string","string", "string", "string"])
            tmp_file = path_obj["uploads"] + "tmp/" + in_json["filename"]
            ext = in_json["filename"].split(".")[-1]
            perform_sitemapping_qc(tmp_file, ext, in_json,mapping_rows,preview_rows)
            field_list = get_field_list(tmp_file, ext)
            file_info = {"filename":in_json["filename"], "fieldlist":field_list}
            out_json = {"mappingrows":mapping_rows, "previewrows":preview_rows,
                    "fileinfo":file_info}
            #Over write tmp_file with additiona qc_comment
            with open(tmp_file, "w") as FW:
                for row in [preview_rows[0]] + preview_rows[2:]: 
                    FW.write("\"%s\"\n" % ("\",\"".join(row)))
        elif in_json["action"] == "save_file":
            in_json["institute"] = in_json["institute"].replace(" ", "").lower()
            tmp_file = path_obj["uploads"] + "tmp/" + in_json["filename"]
            in_json["filename"] = in_json["institute"] + "." + in_json["filename"]
            out_file = path_obj["uploads"] + in_json["filename"]
            cmd = "cp %s %s" % (tmp_file, out_file)
            x = commands.getoutput(cmd)
            meta_file = path_obj["uploads"] + in_json["filename"] + ".json"
            meta_file = meta_file.replace(".csv", "")
            meta_file = meta_file.replace(".tsv", "")
            in_json.pop("action")
            with open(meta_file, "w") as FR:
                FR.write("%s\n" % (json.dumps(in_json, indent=4)))

    except pymongo.errors.ServerSelectionTimeoutError as err:
        out_json = {"taskstatus":0, "errormsg":"Connection to mongodb failed!"}
    except pymongo.errors.OperationFailure as err:
        out_json = {"taskstatus":0, "errormsg":"MongoDB auth failed!"}

    print  json.dumps(out_json, indent=4)



if __name__ == '__main__':
	main()

