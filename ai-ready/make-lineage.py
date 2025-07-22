import os
import sys
import json
import glob
from optparse import OptionParser


import libgly


def get_lineage(go_ns, go_id):

    p_list = []           
    if go_id in go_dict[go_ns]:
        for p in go_dict[go_ns][go_id]["parentlist"]:
            p_list = [p]
            p_list += get_lineage(go_ns, p)

    return p_list


def load_canon2goid(in_file):    
    
    
    canon2goid = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon, go_id = row[0], row[2].replace("_", ":")
        if canon not in canon2goid:
            canon2goid[canon] = []
        if go_id not in canon2goid[canon]:
            canon2goid[canon].append(go_id)

    return canon2goid


def main():

    global go_dict

    ontology_file = "outdir/go-basic.obo"

    is_parent = {}
    obj_list = []
    with open(ontology_file, "r") as FR:
        flag = False
        go_id, go_name, parent_id = "", "", ""
        for line in FR:
            if line[0:6] == "[Term]":
                p_list = []
                flag = True
            if flag == True and line.strip() == "":
                obj_list.append({"id":go_id, "ns":go_ns,"name":go_name, "parentlist":p_list})
                flag = False
            if flag and line[0:3] == "id:":
                go_id = line[3:].strip()
            elif flag and line[0:5] == "name:":
                go_name = line[5:].strip()
            elif flag and line[0:10] == "namespace:":
                go_ns = line[10:].strip()
            elif flag and line[0:5] == "is_a:":
                parent_id = line[5:].strip().split(" ")[0]
                p_list.append(parent_id)
                is_parent[parent_id] = True


    go_dict = {}
    for obj in obj_list:
        go_ns = obj["ns"]
        go_id = obj["id"]
        for parent_id in obj["parentlist"]:
            if go_ns not in go_dict:
                go_dict[go_ns] = {} 
            if go_id not in go_dict[go_ns]:
                go_dict[go_ns][go_id] = {"name":obj["name"], "parentlist":[]}
            if parent_id not in go_dict[go_ns][go_id]["parentlist"]:
                go_dict[go_ns][go_id]["parentlist"].append(parent_id)


    row = ["id","namespace", "lineage"]
    print "\"%s\"" % ("\",\"".join(row))
    for go_ns in go_dict:
        for go_id in go_dict[go_ns]:
            lineage = get_lineage(go_ns, go_id)
            go_name = go_dict[go_ns][go_id]["name"]
            row = [go_id, go_ns, ";".join(lineage)]
            print "\"%s\"" % ("\",\"".join(row))


    
    return


                


if __name__ == '__main__':
        main()
