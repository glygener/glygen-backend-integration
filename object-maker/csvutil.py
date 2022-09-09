import csv
import json
import re

def extractsingleseqsites(df_list, table_ind, action_obj):

    site_aa_list = ["S","T","Y","K","D","E"]

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    peptide_start_field = action_obj["input_fields"]["peptide_start_field"]
    peptide_seq_field = action_obj["input_fields"]["peptide_seq_field"]
    sites_field = action_obj["output_fields"]["sites_field"]

    row_list = []
    for row in data_frame["data"]:
        peptide_seq = row[field_list.index(peptide_seq_field)].upper()
        peptide_start = int(row[field_list.index(peptide_start_field)])
        site_list = []
        for aa in site_aa_list:
            pos_list = [pos.start() for pos in re.finditer(aa, peptide_seq)]
            for pos in pos_list:
                site_list.append("%s%s" % (aa,pos+peptide_start))
        sites = ";".join(site_list)
        newrow = row + [sites]
        row_list.append(newrow)


    new_data_frame = {}
    new_data_frame["fields"] = ""
    new_data_frame["fields"] = field_list + [sites_field]
    new_data_frame["data"] = row_list

    return new_data_frame




def extractpeptideranges(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    sequence_field = action_obj["input_fields"]["sequence_field"]
    peptide_field = action_obj["input_fields"]["peptide_field"]
    start_pos_field = action_obj["output_fields"]["start_pos_field"]
    end_pos_field = action_obj["output_fields"]["end_pos_field"]
    repeat_count_field = action_obj["output_fields"]["repeat_count_field"]

    row_list = []
    for row in data_frame["data"]:
        sequence = row[field_list.index(sequence_field)].upper()
        peptide = row[field_list.index(peptide_field)].upper()
        pos_list = [pos.start() for pos in re.finditer(peptide, sequence)]
        site_count = len(pos_list)
        if site_count == 0:
            newrow = row + ["-1", "-1", "0"]
            row_list.append(newrow)
        else:
            for pos in pos_list:
                start_pos = pos + 1
                end_pos = pos + len(peptide)
                newrow = row + [str(start_pos), str(end_pos), str(site_count)]
                row_list.append(newrow)


    new_data_frame = {}
    new_data_frame["fields"] = field_list + [start_pos_field, end_pos_field, repeat_count_field]
    new_data_frame["data"] = row_list
    
    return new_data_frame


def extractseqsites(df_list, table_ind, action_obj):
   
    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}

    
    row_list = []
    for row in data_frame["data"]:
        pep_start = row[field_list.index(action_obj["start_pos_field"])]
        pep_end = row[field_list.index(action_obj["end_pos_field"])]
        seq = row[field_list.index(action_obj["seq_field"])]
        pep,start_aa, end_aa = "", "", ""
        if pep_start.isdigit() == False or pep_end.isdigit() == False:
            pep = "NA"
            start_aa = "NA"
            end_aa = "NA"
            row_list.append(row + [pep, start_aa, end_aa])
            continue
        if int(pep_start) > int(pep_end) or seq.strip() == "":
            pep = "NA"
            start_aa = "NA"
            end_aa = "NA"
            row_list.append(row + [pep, start_aa, end_aa])
            continue
        if int(pep_start) > len(seq) or int(pep_end) > len(seq):
            pep = "position outside of reference"
            start_aa = "NA"
            end_aa = "NA"
            row_list.append(row + [pep, start_aa, end_aa])
            continue
        pep = seq[int(pep_start)-1:int(pep_end)]
        start_aa, end_aa = pep[0], pep[-1]
        row_list.append(row + [pep, start_aa, end_aa])
    
    data_frame["fields"] += [action_obj["pep_field"]]
    data_frame["fields"] += [action_obj["start_aa_field"]]
    data_frame["fields"] += [action_obj["end_aa_field"]]
    data_frame["data"] = row_list

    return data_frame


def split_col(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
   
    field_idx = field_list.index(action_obj["field"])
    delim = action_obj["delim"]
    new_field_list = field_list[0:field_idx + 1] +  action_obj["newfields"] 
    new_field_list += field_list[field_idx+1:]

    row_list = []
    for row in data_frame["data"]:
        val_list = row[field_idx]
        newrow = row[0:field_idx+1] +  row[field_idx].split(delim) + row[field_idx+1:]
        row_list.append(newrow)
    
    data_frame["data"] = row_list
    data_frame["fields"] = new_field_list

    return data_frame



def add_average_col(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    
    row_grps = {}
    for row in data_frame["data"]:
        key_list = []
        for f in action_obj["anchorfields"]:
            key_list.append(row[field_list.index(f)])
        group_key = " ".join(key_list)
        if group_key not in row_grps:
            row_grps[group_key] = []
        row_grps[group_key].append(row)

    row_list = []
    for group_key in row_grps:
        x_count, x_sum, x_max,x_min = 0.0, 0.0, -100000.0, 1000000.0
        for row in row_grps[group_key]:
            x = float(row[field_list.index(action_obj["field"])])
            x_sum += x
            x_count += 1.0
        x_average = "%.5f" % (x_sum/x_count)
        for row in row_grps[group_key]:
            newrow =  row + [x_average]
            row_list.append(newrow)
    
    data_frame["data"] = row_list
    data_frame["fields"] = field_list + [action_obj["newfield"]]
    return data_frame



    

def add_normalized_col(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    
    row_grps = {}
    for row in data_frame["data"]:
        key_list = []
        for f in action_obj["anchorfields"]:
            key_list.append(row[field_list.index(f)])
        group_key = " ".join(key_list)
        if group_key not in row_grps:
            row_grps[group_key] = []
        row_grps[group_key].append(row)

    row_list = []
    for group_key in row_grps:
        x_sum, x_max,x_min = 0.0, -100000.0, 1000000.0
        for row in row_grps[group_key]:
            x = float(row[field_list.index(action_obj["field"])])
            x_max = x if x > x_max else x_max
            x_min = x if x < x_min else x_min
            x_sum += x

        for row in row_grps[group_key]:
            x = float(row[field_list.index(action_obj["field"])])
            #x_normalized = (x-x_min)/(x_max-x_min)
            x_normalized = "%.5f" % (x/x_sum)
            newrow =  row + [x_normalized]
            row_list.append(newrow)


    data_frame["data"] = row_list
    data_frame["fields"] = field_list + [action_obj["newfield"]]
    return data_frame



def transpose_cols(df_list, table_ind, action_obj):
    
    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    start_col_idx = action_obj["startcolidx"]
    new_field_one = action_obj["newfieldone"]
    new_field_two = action_obj["newfieldtwo"]

    row_list = []
    field_list_new = field_list[0:start_col_idx] + [new_field_one, new_field_two]
    for row in data_frame["data"]:
        for j in range(start_col_idx, len(row)):
            v = row[j]
            newrow = row[0:start_col_idx] + [field_list[j],v]
            row_list.append(newrow)

    data_frame["data"] = row_list
    data_frame["fields"] = field_list_new

    return data_frame


def expand_rows(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
    expansion_field = action_obj["expansionfield"]
    expansion_delim = action_obj["expansiondelim"]
    row_list = []
    for row in data_frame["data"]:
        field_idx = field_list.index(expansion_field)
        val_list = list(set(row[field_idx].strip().split(expansion_delim)))
        
        for val in val_list:
            tmprow = []
            for j in range(0, len(row)):
                tmprow.append(row[j] if j != field_idx else val)
            row_list.append(tmprow)
    
    data_frame["data"] = row_list

    return data_frame


def filter_out_records(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}

    row_list = []
    for row in data_frame["data"]:
        flag_list = []
        for obj in action_obj["conditionlist"]:
            cond_field = obj["field"]
            cond_value = obj["value"]
            data_value = row[field_list.index(cond_field)]
            if obj["operation"] == "in":
                flag_list.append(data_value in cond_value)

        flag_list_unique = list(set(flag_list))
        if len(flag_list) == len(action_obj["conditionlist"]) and flag_list_unique == [True]:
            continue
        row_list.append(row)

    data_frame["data"] = row_list
    return data_frame


def filter_in_records(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}

    row_list = []
    for row in data_frame["data"]:
        flag_list = []
        for obj in action_obj["conditionlist"]:
            cond_field = obj["field"]
            cond_value = obj["value"]
            data_value = row[field_list.index(cond_field)]
            if obj["operation"] == "in":
                flag_list.append(data_value in cond_value)

        flag_list_unique = list(set(flag_list))
        if len(flag_list) == len(action_obj["conditionlist"]) and flag_list_unique == [True]:
            row_list.append(row)

    data_frame["data"] = row_list
    return data_frame



def add_combo_field(df_list, table_ind, action_obj):

    field_list = df_list[table_ind]["fields"] + [action_obj["combofield"]]
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}
   

    for row in data_frame["data"]:
        extra_row = []
        for f in action_obj["fieldlist"]:
            extra_row.append(row[field_list.index(f)].replace("\"", "").strip())
        row += [action_obj["merge_char"].join(extra_row)]

    return data_frame


def add_constant_fields(df_list, table_ind, new_data):

    new_data_field_list = new_data.keys()
    field_list = df_list[table_ind]["fields"] + new_data_field_list
    data_frame = {"fields":field_list, "data":df_list[table_ind]["data"]}

    extra_row = []
    for f in new_data_field_list:
        extra_row.append(new_data[f])
    for row in data_frame["data"]:
        row += extra_row

    return data_frame





def union_tables(df_list,ind_list):

    out_tbl = []
    first_ind = ind_list[0] - 1
    for i in ind_list:
        ind = i - 1
        out_tbl += df_list[ind]["data"]
    
    seen = {}
    data_frame = {"fields":df_list[first_ind]["fields"], "data":[]}
    for row in out_tbl:
        row_str = json.dumps(row)
        if row_str in seen:
            continue
        data_frame["data"].append(row)
        seen[row_str] = True

    return data_frame


def load_large_sheet(sheet_obj, in_file, field_list, separator):

    import sys
    #reload(sys)
    #sys.setdefaultencoding('utf-8')
    import io

    sheet_obj["fields"] = []
    sheet_obj["data"] = []
    seen = {}
    field_ind_list = []
    
    #with io.open(in_file, "r") as FR:
    with io.open(in_file, "r", encoding="utf-8-sig",errors="ignore") as FR:
        rowcount = 0
        ncols = 0
        for line in FR:
            row = line.strip().split(separator)
            #row[0] = row[0].split("\"")[1]
            row[0] = row[0].replace("\"", "")
            row[-1] = row[-1].replace("\"", "")
            rowcount += 1
            if rowcount == 1:
                for j in range(0, len(row)):
                    row[j] = row[j].strip().replace("\"", "")
                bad_fields = []
                for f in field_list:
                    if f not in row:
                        bad_fields.append(f)
                if bad_fields != []:
                    print ("input file:",  in_file)
                    print ("fields in file:", row)
                    print ("fields for extraction:", field_list)
                    print ("non matching fields: ", bad_fields)
                    sys.exit()

                
                #capture number of columns here
                ncols = len(row)
                field_list = row if field_list == [] else field_list
                for f in field_list:
                    field_index = row.index(f)
                    if field_index != -1:
                        sheet_obj["fields"].append(f)
                        field_ind_list.append(field_index)
            else:
                #make sure every row has ncols columns
                if len(row) != ncols:
                    print ("bad row %s" % (rowcount))
                    print (row)
                    sys.exit()
                new_row = []
                for j in field_ind_list:
                    new_row.append(row[j].strip())
                if json.dumps(new_row) not in seen:
                    sheet_obj["data"].append(new_row)
                    seen[json.dumps(new_row)] = True
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
            if json.dumps(row) in seen:
                continue
            seen[json.dumps(row)] = True
            row_count += 1
            if row_count == 1:
                f_list = row
                anchor_field = row[0] if anchor_field == "" else anchor_field
                for j in range(0, len(row)):
                    if row[j] == anchor_field:
                        continue
                    sheet_obj["fields"].append(row[j].strip().replace("\"", ""))
            else:
                new_row = []
                for j in range(0, len(row)):
                    if f_list[j] == anchor_field:
                        continue
                    new_row.append(row[j].strip().replace("\"", "`"))
                main_id = row[f_list.index(anchor_field)]
                if main_id not in sheet_obj["data"]:
                    sheet_obj["data"][main_id] = [] 
                sheet_obj["data"][main_id].append(new_row)
    return


def load_sheet_from_fasta(sheet_obj, in_file, field_list):

    from Bio import SeqIO

    sheet_obj["fields"] = field_list
    sheet_obj["data"] = []
    for record in SeqIO.parse(in_file, "fasta"):
        seq_id = record.id
        if record.id.find("|") != -1:
            seq_id = record.id.split("|")[1]
        seq_desc = " ".join(record.description.split(" ")[1:])
        row = [seq_id, seq_desc, str(record.seq.upper())]
        sheet_obj["data"].append(row)

    return




def load_sheet_from_json(sheet_obj, in_file, field_list):

    sheet_obj["fields"] = field_list
    sheet_obj["data"] = []
    doc = json.loads(open(in_file, "r").read())
    for obj in doc:
        val_dict = {}
        for prop_lineage in field_list:
            val = obj
            for prop in prop_lineage.split("."):
                if prop in val:
                    val = val[prop]
                else:
                    val = ""
            val_dict[prop_lineage] = val
        row = []
        for f in field_list:
            val = str(val_dict[f]) if f in val_dict else ""
            row.append(val)
        sheet_obj["data"].append(row)
    
    return

def get_sheet_stats(in_file, separator):

    seen_row, seen_id = {}, {}
    row_count, id_count, field_count = 0, 0, 0
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=",", quotechar='\"')
        idx = 0
        for row in csv_grid:
            if json.dumps(row) in seen_row:
                continue
            seen_row[json.dumps(row)] = True
            if idx == 0:
                field_count = len(row)
            else:
                row_count += 1
                seen_id[row[0]] = True
            idx += 1
    id_count = len(list(seen_id.keys()))

    return field_count, row_count, id_count


def load_sheet(sheet_obj, in_file, field_list, separator):

    seen = {}
    sheet_obj["fields"] = []
    sheet_obj["data"] = []
    field_ind_list = []
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=",", quotechar='\"')
        row_count = 0
        ncols = 0
        for row in csv_grid:
            if json.dumps(row) in seen:
                continue
            seen[json.dumps(row)] = True
            row_count += 1
            for j in range(0, len(row)):
                row[j] = row[j].replace("\"", "`")
            if row_count == 1:
                ncols = len(row)
                for j in range(0, len(row)):
                    f = row[j].strip()
                    if field_list != [] and f in field_list:
                        field_ind_list.append(j)
                        sheet_obj["fields"].append(f)
                    else:
                        field_ind_list.append(j)
                        sheet_obj["fields"].append(f)
            else:
                #make sure every row has ncols columns
                if len(row) != ncols:
                    continue
                new_row = []
                for j in field_ind_list:
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



def join_tables(df_one, df_two, anchor_fields):
    #It is assumed that the first column in both tables are the anchor columns
    
    df_three = {"fields":[], "data":[]}
    df_three["fields"] = []
    for f in df_one["fields"]:
        df_three["fields"].append(f)

    for f in df_two["fields"]:
        if f != anchor_fields[1]:
            df_three["fields"].append(f)

    row_dict_two = {}
    for row in df_two["data"]:
        anchor_id = row[df_two["fields"].index(anchor_fields[1])]
        if anchor_id not in row_dict_two:
            row_dict_two[anchor_id] = []
        newrow = []
        for j in range(0, len(df_two["fields"])):
            if df_two["fields"][j] != anchor_fields[1]:
                newrow.append(row[j])
        row_dict_two[anchor_id].append(newrow)

    row_two_empty = []
    for j in range(1, len(df_two["fields"])):
        row_two_empty.append("")


    for row_one in df_one["data"]:
        anchor_id = row_one[df_one["fields"].index(anchor_fields[0])]
        if anchor_id in row_dict_two:
            for row_two in row_dict_two[anchor_id]:
                df_three["data"].append(row_one + row_two)
        else:
            df_three["data"].append(row_one + row_two_empty)


    return df_three







