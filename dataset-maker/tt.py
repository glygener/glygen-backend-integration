import sys
import csv
import glob
import gzip


import libgly



def main():

    in_file = "downloads/glytoucan/current/export/glycoctxml.zip.00"
    with open(in_file, 'r') as f:
         print(f.readlines())
    exit()
 
    bin_list = [{"s":0, "e":0}]
    for i in range(0, 100):
        o = {"s":i*100 + 1, "e":(i+1)*100}
        bin_list.append(o)

    count_dict = {}
    with open("tmp/junk.txt", "r") as FR:
        for line in FR:
            n = int(line.strip())
            for i in range(0,len(bin_list)):
                s, e = bin_list[i]["s"], bin_list[i]["e"]
                if i not in count_dict:
                    count_dict[i] = 0
                if n == 0:
                    count_dict[i] += 1
                elif n >= s and n <= e:
                    count_dict[i] += 1

    for i in range(0,len(bin_list)):
        print count_dict[i], "MissingScore in [%s-%s]" % (bin_list[i]["s"], bin_list[i]["e"])


if __name__ == '__main__':
    main()




