import os
import csv
import json
import time
import subprocess



def log_file_usage(used_file, ds, flag):

    pid = os.getpid()
    log_file = "usage/file_usage.%s.log" % (pid)

    line = used_file + "\n" if used_file != "" else ""
    FW = open(log_file, "a") if flag == "append" else open(log_file, "w")
    FW.write("%s" % (line))
    FW.close()

    return



def log_usage(log_file, mem_before, start_time):

    mem_after = process_memory()
    elapsed_time = time.time() - start_time
    f = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    #with open(log_file, "w") as FW:
    #    FW.write("%s, elapsed time = %s, memory usage = %s\n" % (log_file,f,mem_after-mem_before))

    return


def process_memory():
    #process = psutil.Process(os.getpid())
    #mem_info = process.memory_info()
    #return mem_info.rss
    return {}



def load_species_info(species_obj, in_file):

    data_frame = {}
    load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        obj = {}
        for f in f_list:
            obj[f] = row[f_list.index(f)]
        tax_id = obj["tax_id"] 
        short_name = obj["short_name"]
        if tax_id not in species_obj:
            species_obj[tax_id] = {}
            species_obj[short_name] = {}
        for f in obj:
            species_obj[tax_id][f] = int(obj[f]) if f == "tax_id" else obj[f] 
            species_obj[short_name][f] = int(obj[f]) if f == "tax_id" else obj[f]

    return


def load_sheet_as_dict(sheet_obj, in_file, separator, anchor_field):


    seen = {}
    
    if "fields" not in sheet_obj:
        sheet_obj["fields"] = []
    if "data" not in sheet_obj:
        sheet_obj["data"] = {}

    f_list = []
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='\"')
        row_count = 0
        for row in csv_grid:
            if row == []:
                continue
            if json.dumps(row) in seen:
                continue
            seen[json.dumps(row)] = True
            row_count += 1
            if row_count == 1:
                f_list = row
                for j in range(0, len(row)):
                    if row[j] == anchor_field:
                        continue
                    sheet_obj["fields"].append(row[j].strip().replace("\"", ""))
            else:
                if len(f_list) != len(row):
                    continue
                new_row = []
                for j in range(0, len(row)):
                    if f_list[j] == anchor_field:
                        continue
                    new_row.append(row[j].strip())
                main_id = row[f_list.index(anchor_field)]
                if main_id not in sheet_obj["data"]:
                    sheet_obj["data"][main_id] = [] 
                sheet_obj["data"][main_id].append(new_row)
    return



def load_sheet(sheet_obj, in_file, separator):


    cmd = "readlink -f " + in_file
    x = subprocess.getoutput(cmd)
    log_file_usage(x, "", "append")


    seen = {}
    sheet_obj["fields"] = []
    sheet_obj["data"] = []
    #with io.open(in_file, "r", encoding="utf-8-sig",errors="ignore") as FR:
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='"')
        row_count = 0
        for row in csv_grid:
            if json.dumps(row) in seen:
                continue
            seen[json.dumps(row)] = True
            row_count += 1
            if row_count == 1:
                for j in range(0, len(row)):
                    sheet_obj["fields"].append(row[j].strip().replace("\"", ""))
            else:
                new_row = []
                for j in range(0, len(row)):
                    new_row.append(row[j].strip())
                sheet_obj["data"].append(new_row)
    return


def load_workbook(workbook_obj, fileset_objlist, separator):

    for obj in fileset_objlist:
        for file_name in obj["filenamelist"]:
            in_file = obj["dir"] + file_name
            workbook_obj["sheets"][file_name] = {}
            load_sheet(workbook_obj["sheets"][file_name], in_file, ",")

    return

def left_join_tables(tbl_a, tbl_b):

    map_dict = {}
    for i in range(0, len(tbl_b)):
        if tbl_b[i][0] not in map_dict:
            map_dict[tbl_b[i][0]] = []
        map_dict[tbl_b[i][0]].append(i)

    empty_row = []
    for j in range(1, len(tbl_b[0])):
        empty_row.append("")

    tbl_c = [tbl_a[0] + tbl_b[0][1:]]
    for row in tbl_a[1:]:
        main_id_a = row[0]
        if main_id_a in map_dict:
            for i in map_dict[main_id_a]:
                tbl_c.append(row + tbl_b[i][1:])
        else:
            tbl_c.append(row + empty_row)

    return tbl_c


def get_doi_citation(doi, in_file):

    row = []
    if os.path.isfile(in_file) == True:
        data_frame = {}
        load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            if row[f_list.index("doi_id")] == doi:
                newrow = [
                    row[f_list.index("title")],
                    row[f_list.index("journal_name")],
                    row[f_list.index("publication_date")],
                    row[f_list.index("authors")]
                ]
                return newrow

    return []



def get_citation(pmid):
    
    
    in_file_one = "downloads/ncbi/medline_txt/pmid.%s.txt" % (pmid)
    in_file_two = "downloads/ncbi/medline_json/pmid.%s.json" % (pmid)
    obj = {}
    if os.path.isfile(in_file_one) == False and os.path.isfile(in_file_two) == False:
        return {"row":[], "flaglist":["not-downloaded"]}
    elif os.path.isfile(in_file_one) == True:
        line_list = open(in_file_one, "r").read().split("\n")
        lcount = 0
        prev_key = ""
        for line in line_list:
            lcount += 1
            if lcount > 3:
                key = line[0:4].strip()
                val = line[5:].strip()
                if key not in obj:
                    obj[key] = []
                if key == "":
                    obj[prev_key].append(val)
                else:
                    obj[key].append(val)
                    prev_key = key
    elif os.path.isfile(in_file_two) == True:
        obj = json.loads(open(in_file_two, "r").read())

    if "AU" not in obj and "IR" in obj:
        obj["AU"] = obj["IR"]
    if "AU" not in obj:
        obj["AU"] = ["author information missing"]
    for tag in ["TI", "DP", "JT", "BTI"]:
        if tag not in obj:
            obj[tag] = ["missing tag=%s" % (tag)]

    row, flag_list = [], []
    #for tag in ["TI", "DP"]:
    #    if tag not in obj:
    #        flag_list.append("missing-%s" % (tag))
    #if "JT" not in obj and "BTI" not in obj:
    #    flag_list.append("missing-%s-or-%s" % ("JT", "BTI"))
        
    if flag_list == []:
        if type(obj["TI"]) is list:
            title = " ".join(obj["TI"]).replace("\"", "`")
            journal = " ".join(obj["JT"]) if "JT" in obj else " ".join(obj["BTI"])
            pubdate = " ".join(obj["DP"])
            authors = ", ".join(obj["AU"])
            row = [title, journal, pubdate, authors]
        else:
            title = obj["TI"].replace("\"", "`")
            journal = obj["JT"] if "JT" in obj else obj["BTI"]
            pubdate = obj["DP"]
            authors = obj["AU"]
            row = [title, journal, pubdate, authors]

    return {"row":row, "flaglist":flag_list}

