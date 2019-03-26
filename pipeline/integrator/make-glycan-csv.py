import os,sys
import json
import glob
from optparse import OptionParser

__version__="1.0"
__status__ = "Dev"

###############################
def load_mapping(map_dict, in_file, map_name,header_list):

    map_dict[map_name] = {}
    with open(in_file, "r") as FR:
        row_count = 0
        for line in FR:
            row_count += 1
            if row_count == 1:
                header_list.append(line.strip().split("\t")[1])
                continue
            ids = line.strip().split("\t")
            if len(ids)==2:
                id1=ids[0]
                id2=ids[1]
            else:
                id1=ids[0]
                id2='-'
            if id1 not in map_dict[map_name]:
                map_dict[map_name][id1] = []
            map_dict[map_name][id1].append(id2)
    return
###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-o","--organismType",action="store",dest="organism_type",help="The type of organism (human, mouse)")
    parser.add_option("-L","--List",action="store",dest="List",help='List'+" '"+'["pbch","taxa"]'+"'")

    (options,args) = parser.parse_args()
    for file in ([options.organism_type, options.List]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    organism = options.organism_type
    List = json.loads(options.List)
    config_json = json.loads(open("conf/config-1.json", "r").read())
    in_dir = config_json["pathinfo"]["glytoucan_input"]
    pattern = in_dir + "*"
    input_file_list = [
        in_dir+"human_gtc6.txt"
        ,in_dir+"mouse_gtc6.txt"
    ]

    file_lists = glob.glob(pattern)
    temp = []
    for element in List:
        for file_list in file_lists:
            txt_file_name = file_list.split('/')[-1]
            if element in txt_file_name:
                temp.append(element)
                input_file_list.append(file_list)
                break
    temp = '_'.join(temp)
    out_file = "/data/share/"+ "/"+organism+"_glycan_"+temp+".csv" #config_json["pathinfo"]["glytoucan_output"]+"/"+organism+"_glycan_"+temp+".csv"

    seen = {"gtid":{}}
    if organism=="human":
        input_file = input_file_list[0]
    else:
        input_file = input_file_list[1]
    with open(input_file, "r") as FR:
        row_count = 0
        for line in FR:
            row_count += 1
            if row_count == 1:
                continue
            gt_id = line.strip().split("\t")[0]
            seen["gtid"][gt_id] = True

    header_list = ["glytoucan_acc"]
    for idx, val in enumerate(List):
        load_mapping(seen, input_file_list[idx+2], val, header_list)

    fw = open(out_file, "w")
    fw.write("%s\n" % (",".join(header_list)))
    for gt_id in seen["gtid"]:
        temp = [gt_id]
        for element in List:
            if gt_id in seen[element]:
                if '-' in seen[element][gt_id]:
                    List1 = [i for i in set(seen[element][gt_id])]
                    if '-' in List1:
                        temp.append("")
                    else:
                        temp.append("|".join(List1))
                else:
                    temp.append("|".join(seen[element][gt_id]))
            else:
                temp.append("")
        fw.write("%s\n" % (",".join(temp)))
    fw.close()

if __name__ == '__main__':
    main()
