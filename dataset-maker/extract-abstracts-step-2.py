import os
import sys
import json
import glob


def main():

    DEBUG = False
    #DEBUG = True


    xml_dir = "downloads/ncbi/medline_xml/"
    xml_file_list = sorted(glob.glob(xml_dir + "/*.xml.gz"))
    n = len(xml_file_list)
    batch_size = int(len(xml_file_list)/10) if len(xml_file_list) > 10 else len(xml_file_list)
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
        cmd = "nohup python3 batch-extract-abstracts.py -s %s -e %s  & " % (o["s"], o["e"])
        os.system(cmd)
        #print (cmd)



if __name__ == '__main__':
    main()


