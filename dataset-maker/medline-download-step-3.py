import os
import sys
import gzip
import json
from optparse import OptionParser
import glob
import pubmed_parser as pp






###############################
def main():

    doc = json.loads(open("conf/medline.json", "r").read())

    category = "updatefiles"
    updatefiles_start = doc[category]["start"]
   

    xml_folder = "downloads/ncbi/medline_xml/"
    out_folder = "downloads/ncbi/medline_list/"
    tmp_list = glob.glob("downloads/ncbi/medline_xml/*.gz")    
    gz_xml_file_list = []
    for in_file in tmp_list:
        idx = in_file.split("pubmed25n")[-1].replace(".xml.gz", "")
        if int(idx) >= updatefiles_start:
            gz_xml_file_list.append(in_file)
    gz_xml_file_list = sorted(gz_xml_file_list)
    #gz_xml_file_list = sorted(tmp_list)
  
    with open("logs/medline_xml_file_list.json", "w") as FW:
        FW.write("%s\n" % (json.dumps(gz_xml_file_list)))

 
    n = len(gz_xml_file_list)
    batch_size = int(len(gz_xml_file_list)/10) if len(gz_xml_file_list) > 10 else len(gz_xml_file_list)

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
        cmd = "nohup python3 extract-pmid-list.py -s %s -e %s  & " % (o["s"], o["e"])
        os.system(cmd)
        #print (cmd)





if __name__ == '__main__':
    main()

