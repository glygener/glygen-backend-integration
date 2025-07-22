import os
import sys
import gzip
import json
from optparse import OptionParser
import glob
import pubmed_parser as pp






###############################
def main():

    file_list = glob.glob("alignments/homologset/*.aln")
    file_list += glob.glob("alignments/isoformset/*/*.aln")
    file_list = sorted(file_list)

    n = len(file_list)
    batch_size = int(len(file_list)/10) if len(file_list) > 10 else len(file_list)
    range_list = [] 
    for i in range(0, 12):
        s = i*batch_size + 1
        e = s + batch_size
        if e >= n:
            e = n
            range_list.append({"s":s, "e":e})
            break
        range_list.append({"s":s, "e":e}) 

    for o in range_list:
        cmd = "nohup python3 alignmentdb-batch.py -s %s -e %s  & " % (o["s"], o["e"])
        os.system(cmd)
        #print (cmd)





if __name__ == '__main__':
    main()

