import os,sys
import json
import csv
import glob

__version__="1.0"
__status__ = "Dev"



###############################
def main():

    idx = 1
    row_list = [] 
    in_file = "generated/datamodel/datamodel.csv"
    with open(in_file, "r") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            if line[0] == "#":
                continue
            if line.strip() == "":
                continue
            parts = line.strip().split(",")
            prd, dmn, rng, cat = parts[0], parts[1], parts[2], parts[3]
            if rng in ["xsd:string", "xsd:int"]:
                rng = "rdfs:Literal"
            row_list.append([prd, dmn, rng, cat])
    

    out_file = "generated/datamodel/datamodel_decoupled.csv"
    with open (out_file, "w") as FW:
        for row in row_list:
            FW.write("%s\n" % (",".join(row)))




if __name__ == '__main__':
	main()
