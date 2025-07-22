import os,sys
import string
import glob
import json
import subprocess



###############################
def main():



    # Make sure to update "conf/medline.json"
    doc = json.loads(open("conf/medline.json", "r").read())

    xml_folder = "downloads/ncbi/medline_xml/"
    zeros = "0000000"
    #cat_list = ["baseline", "updatefiles"]
    cat_list = ["updatefiles"]
    for category in cat_list:
        start = doc[category]["start"]
        end = doc[category]["end"]
        year = doc[category]["year"]
        ftp_url = doc[category]["url"]
        for idx in range(start, end + 1):
            padlen = 4 - len(str(idx))
            idx = zeros[:padlen] + str(idx)
            file_name = "pubmed%sn%s.xml.gz" % (year, idx)
            out_file_one = xml_folder + file_name
            #if os.path.isfile(out_file_one) == True:
            #    continue
            cmd = "wget %s%s -O %s" % (ftp_url, file_name, out_file_one)
            x = subprocess.getoutput(cmd)
            #print (cmd)
        



if __name__ == '__main__':
	main()

