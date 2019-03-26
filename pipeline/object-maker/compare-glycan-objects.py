import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


###############################
def main():

    in_file_one = "/data/projects/glygen/releases/data/v-1.0.3/jsondb/glycandb.json"
    in_file_two = "/data/projects/glygen/jsondb/glycandb.json"
                
    doc_list_one = json.loads(open(in_file_one, "r").read()) 
    doc_list_two = json.loads(open(in_file_two, "r").read())


    for ac in doc_list_one:
        for o in doc_list_one[ac]["glycoprotein"]:
            o["position"] = int(o["position"])
    
    #for ac in doc_list_two:
    #    print "A|", ac, ac in doc_list_one
    #for ac in doc_list_one:
    #    print "B|", ac, ac in doc_list_two

    for ac in doc_list_two:
        if ac in doc_list_one:
            if doc_list_one[ac] != doc_list_two[ac]:
                #for k in doc_list_two[ac]:
                #    print k, doc_list_one[ac][k] == doc_list_two[ac][k]
                #continue

                for k in ["enzyme", "crossref"]:
                    if doc_list_one[ac][k] != doc_list_two[ac][k]:
                        #print json.dumps(doc_list_one[ac][k], indent=4)
                        #print json.dumps(doc_list_two[ac][k], indent=4)

                        
                        key_list = doc_list_one[ac][k][0].keys()
                        for key in key_list:
                            if key == "evidence":
                                continue
                            set_one, set_two = [], []
                            for j in xrange(0, len(doc_list_one[ac][k])):
                                set_one.append(doc_list_one[ac][k][j][key])
                                if j < len(doc_list_two[ac][k]):
                                    set_two.append(doc_list_two[ac][k][j][key])
                            s_one = sorted(set(set_one))
                            s_two = sorted(set(set_two))
                            print ac, k, key, s_one == s_two, len(s_one), len(s_two)
                            #if sorted(set(set_one)) !=  sorted(set(set_two)):
                            #    print ac, key, sorted(set(set_one)), sorted(set(set_two))








if __name__ == '__main__':
	main()

