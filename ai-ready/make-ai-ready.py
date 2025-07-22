import os
import sys
import json
import glob
from optparse import OptionParser


import libgly

def load_goid2name():

    goid2name = {}
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
                flag = False
            if flag and line[0:3] == "id:":
                go_id = line[3:].strip()
            elif flag and line[0:5] == "name:":
                go_name = line[5:].strip()
                goid2name[go_id] = go_name
    return goid2name


def load_data_matrix(in_file, transpose):

    seen = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        row_id = row[0] if transpose == False else row[1]
        col_id = row[1] if transpose == False else row[0]
        if row_id not in seen:
            seen[row_id] = {}
        seen[row_id][col_id] = 1 if len(f_list) == 2 else float(row[2])

    return seen


def load_goid2lineage(in_file, namespace):

    goid2lineage = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        go_id, go_ns, lineage = row[0], row[1], row[2].split(";")
        if namespace != "" and namespace != go_ns:
            continue
        if go_id not in goid2lineage:
            goid2lineage[go_id] = {}
        goid2lineage[go_id][go_ns] = lineage

    return goid2lineage


def load_canon2goid(species, go_ns, depth):

    lineage_file = "outdir/go_lineage.csv"
    ann_file = "unreviewed/%s_protein_go_annotation.csv" % (species)
     
    goid2lineage = load_goid2lineage(lineage_file, go_ns)

    tmp_dict = {}
    data_frame = {}
    delim = "," if ann_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, ann_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        key, go_id = row[0], row[2].replace("_", ":")
        if go_id not in goid2lineage:
            continue
        for go_ns in goid2lineage[go_id]:
            lineage = goid2lineage[go_id][go_ns]
            if len(lineage) > depth:
                last_p = lineage[-depth]
                value = last_p + "|" + go_ns
                if key not in tmp_dict:
                    tmp_dict[key] = []
                if value not in tmp_dict[key]:
                    tmp_dict[key].append(value)
    return tmp_dict




def load_simplemap(in_file, k_field, v_field):

    tmp_dict = {}
    data_frame = {}
    delim = "," if in_file.split(".")[-1] == "csv" else "\t"
    libgly.load_sheet(data_frame, in_file, delim)
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        key = row[f_list.index(k_field)].strip()
        value = row[f_list.index(v_field)].strip()
        if key == "" or value == "":
            continue
        if key not in tmp_dict:
            tmp_dict[key] = []
        if value not in tmp_dict[key]:
            tmp_dict[key].append(value)

    return tmp_dict




def collapse_matrix_cols(in_matrix, map_dict):
    count_dict = {}
    for row_id in in_matrix:
        for col_id in in_matrix[row_id]:
            val = in_matrix[row_id][col_id]
            if col_id not in map_dict:
                continue
            if map_dict[col_id] == []:
                map_dict[col_id] = ["other"]

            for new_col_id in map_dict[col_id]:
                if row_id not in count_dict:
                    count_dict[row_id] = {}
                if new_col_id not in count_dict[row_id]:
                    count_dict[row_id][new_col_id] = 0
                count_dict[row_id][new_col_id] += val
    
    return count_dict




def collapse_matrix_rows(in_matrix, map_dict):
    
    col_id_list = []
    for row_id in in_matrix:
        for col_id in in_matrix[row_id]:
            if col_id not in col_id_list:
                col_id_list.append(col_id)

    count_dict = {}
    for row_id in in_matrix:
        if row_id not in map_dict:
            continue
        if map_dict[row_id] == []:
            map_dict[row_id] = ["other"]

        for col_id in col_id_list:
            val = in_matrix[row_id][col_id] if col_id in in_matrix[row_id] else 0
            for new_row_id in map_dict[row_id]:
                if new_row_id not in count_dict:
                    count_dict[new_row_id] = {}
                if col_id not in count_dict[new_row_id]:
                    count_dict[new_row_id][col_id] = 0
                count_dict[new_row_id][col_id] += val
    
    return count_dict




def normalize_matrix(in_matrix):
    
    sum_dict = {}
    for row_id in in_matrix:
        for col_id in in_matrix[row_id]:
            if col_id not in sum_dict:
                sum_dict[col_id] = 0.0
            sum_dict[col_id] += 1.0

    for row_id in in_matrix:
        for col_id in in_matrix[row_id]:
            in_matrix[row_id][col_id] = round(float(in_matrix[row_id][col_id])/float(sum_dict[col_id]), 8)

    return sorted(list(sum_dict.keys()))

def get_col_id_list(in_matrix):

    sum_dict = {}
    for row_id in in_matrix:
        for col_id in in_matrix[row_id]:
            if col_id not in sum_dict:
                sum_dict[col_id] = 0.0
            sum_dict[col_id] += 1.0

    return sorted(list(sum_dict.keys()))

def get_matrix_dim(in_matrix):

    seen_row, seen_col = {}, {}
    for row_id in in_matrix:
        seen_row[row_id] =  True
        for col_id in in_matrix[row_id]:
            seen_col[col_id] = True

    return len(list(seen_row.keys())), len(list(seen_col.keys()))


def dump_enrichment(data_matrix, row_map):


    goid2name = load_goid2name()
    seen_one, seen_two = {}, {}
    canon_count = 0
    for canon in row_map:
        for go_id in row_map[canon]:
            if go_id not in seen_one:
                seen_one[go_id] = {}
                seen_two[go_id] = {}
            if canon in data_matrix:
                seen_two[go_id][canon] = True
            seen_one[go_id][canon] = True
        canon_count += 1




    for go_id in seen_one:
        n_one = len(list(seen_one[go_id].keys()))
        n_two = len(list(seen_two[go_id].keys()))
        go_id = go_id.split("|")[0]
        name = goid2name[go_id] if go_id in goid2name else "unknown"
        r = round(float(n_two)/float(n_one), 4)
        print ("%s (%s/%s) | %s (%s) " % (r,n_two, n_one,name, go_id))

    return


def load_gtc2species():

    in_file = "outdir/species2glycan.csv"
    mat = load_data_matrix(in_file, True)
    tmp_dict = {}
    idx = 1
    seen = {}
    count_dict = {}
    for gtc in mat:
        sp_list = sorted(list(mat[gtc].keys()))
        if len(sp_list) > 1:
            continue
        sp = ",".join(sp_list)
        if sp not in seen:
            seen[sp] = idx
            idx += 1
        tmp_dict[gtc] = {"idx":seen[sp], "sp":sp}
        if sp not in count_dict:
            count_dict[sp] = 0
        count_dict[sp] += 1

    final_dict = {}
    for gtc in tmp_dict:
        idx, sp = tmp_dict[gtc]["idx"], tmp_dict[gtc]["sp"]
        if count_dict[sp] > 100:
            final_dict[gtc] = {"idx":idx, "sp":sp}
    
    return final_dict


def main():


    #species = "human"
    species = "mouse"
    #species = "rat"
    
    #dataset = "glycan2protein"
    #dataset = "motif2glycan"
    #dataset = "enzyme2glycan"
    #dataset = "species2glycan"
    dataset = "glycan2compositionadv"

    row_map, col_map = {}, {}

    #in_file = "outdir/%s_%s.csv" % (species, dataset)
    in_file = "outdir/%s.csv" % (dataset)
    #transpose = True
    transpose = False
    data_matrix = load_data_matrix(in_file, transpose)


    gtc2species = load_gtc2species()


    # "molecular_function", "biological_process", "cellular_component"
    #go_ns, depth = "cellular_component", 3
    #col_map = load_canon2goid(species, go_ns, depth)
    #dump_enrichment(data_matrix, row_map)
    #exit()

    #in_file = "unreviewed/glycan_classification.csv"
    #k_field, v_field = "glytoucan_ac", "glycan_subtype"
    #row_map = load_simplemap(in_file, k_field, v_field)
    #in_file = "unreviewed/glycan_enzyme.csv"
    #k_field, v_field = "glytoucan_ac", "uniprotkb_canonical_ac"
    #row_map = load_simplemap(in_file, k_field, v_field)
    
    #in_file = "unreviewed/glycan_motif.csv"
    #k_field, v_field = "glytoucan_ac", "motif_ac"
    #col_map = load_simplemap(in_file, k_field, v_field)
    

    #in_file = "unreviewed/glycan_enzyme.csv"
    #k_field, v_field = "glytoucan_ac", "uniprotkb_canonical_ac"
    #row_map = load_simplemap(in_file, k_field, v_field)


    if row_map != {}:
        data_matrix = collapse_matrix_rows(data_matrix, row_map)
    if col_map != {}:
        data_matrix = collapse_matrix_cols(data_matrix, col_map)
    

    #col_id_list = normalize_matrix(data_matrix)
    col_id_list = get_col_id_list(data_matrix)

    row = [""] + col_id_list + ["class"]
    print ("\"%s\"" % ("\",\"".join(row)))
    for row_id in data_matrix:
        if row_id not in gtc2species:
            continue
        row = [row_id]
        sp_class = str(gtc2species[row_id]["idx"])
        for col_id in col_id_list:
            r = str(data_matrix[row_id][col_id]) if col_id in col_id in data_matrix[row_id] else "0.0"
            row.append(r)
        print ("\"%s\"" % ("\",\"".join(row + [sp_class])))

            
    return


                


if __name__ == '__main__':
        main()
