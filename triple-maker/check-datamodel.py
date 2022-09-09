import os,sys
import json
import csv
import glob


__version__="1.0"
__status__ = "Dev"


def load_ann(in_file):

    ann = {}
    with open(in_file, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_frame:
            row_count += 1
            if row_count > 1:
                cls = row[0]
                if cls not in ann:
                    ann[cls] = {
                        "label":row[1], "comment":row[2], "definedby":row[3], "parent":row[4], "seealso":row[5]
                    }
    return ann


###############################
def main():

    
    config_obj = json.loads(open("../conf/config.json", "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
    ns_map =  config_obj["nsmap"]
    uri_map = config_obj["urimap"]

    seen = {"prd_in_data":{}, "prd_in_model":{}}
    
    in_file = "generated/datamodel/datamodel.csv"
    with open(in_file, "r") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            if lcount <= 1 or line[0] == "#":
                continue
            if line.strip() == "":
                continue
            parts = line.strip().split(",")
            prd_one = parts[0]
            ns, prd_two = parts[0].split(":")
            prd = prd_one
            if ns in ns_map:
                prd = ns_map[ns] + prd_two
            prd = prd.replace("www.", "")
            seen["prd_in_model"][prd] = True


    file_list = glob.glob("outdir/glygen/*.nt")
    for in_file in file_list:
        with open(in_file, "r") as FR:
            for line in FR:
                parts = line.strip().split(" ")
                prd = parts[1].strip()[1:-1]
                prd = prd.replace("www.", "")
                seen["prd_in_data"][prd] = True



    for prd in seen["prd_in_data"]:
        if prd in seen["prd_in_model"]:
            print "DM,%s" %(prd)
        else:
            print "D0,%s" %(prd)

    for prd in seen["prd_in_model"]:
        if prd not in seen["prd_in_data"]:
            print "M0,%s" %(prd)



    sys.exit()



if __name__ == '__main__':
	main()
