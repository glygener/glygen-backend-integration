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


def compare_objects(main_id, obj_one, obj_two, ind_name):

    if obj_one != obj_two:
        if type(obj_one) is list:
            if len(obj_one) != len(obj_two):
                print "lendiff|%s|%s|%s|%s" % (main_id, ind_name, len(obj_one), len(obj_two))
            else:
                for i in xrange(0, len(obj_one)):
                    new_ind_name = ind_name + "_%s" % (i) 
                    if type(obj_one[i]) is list:
                        obj_one[i] = sorted(obj_one[i])
                        obj_two[i] = sorted(obj_two[i])
                    compare_objects(main_id, obj_one[i], obj_two[i], new_ind_name)
        elif type(obj_one) is dict:
            for k in obj_one:
                new_ind_name = ind_name + "_%s" % (k) 
                if k in obj_two:
                    compare_objects(main_id, obj_one[k], obj_two[k], new_ind_name)
        else:
            print "valdiff|%s|%s|%s|%s" % (main_id, ind_name, obj_one, obj_two)




###############################
def main():

    in_file_one = "/data/projects/glygen/releases/data/v-1.0.3/jsondb/proteindb.json"
    in_file_two = "/data/projects/glygen/jsondb/proteindb.json"
                
    doc_list_one = json.loads(open(in_file_one, "r").read()) 
    doc_list_two = json.loads(open(in_file_two, "r").read())


    k = "pathway"
    for ac in doc_list_one:
        if ac in doc_list_two:
            if type(doc_list_one[ac][k]) is list:
                doc_list_one[ac][k] = sorted(doc_list_one[ac][k])
                doc_list_two[ac][k] = sorted(doc_list_two[ac][k])
            compare_objects(ac, doc_list_one[ac][k], doc_list_two[ac][k], k)
    sys.exit()




    #for ac in doc_list_two:
    #    if ac not in doc_list_one:
    #        print "new,%s,%s" % (ac, doc_list_two[ac]["species"][0]["name"])

    #for ac in doc_list_one:
    #    if ac not in doc_list_two:
    #        print "retired,%s,%s" % (ac, doc_list_one[ac]["species"][0]["name"])
    #sys.exit()


    for ac in doc_list_two:
        if ac in doc_list_one:
            if doc_list_one[ac] != doc_list_two[ac]:
                for k in doc_list_two[ac]:
                    if doc_list_one[ac][k] != doc_list_two[ac][k]:
                        print "%s,%s,%s" % (ac, k, doc_list_two[ac]["species"][0]["name"])
                continue

                for k in ["mass"]:
                    if doc_list_one[ac][k] != doc_list_two[ac][k]:
                        #print json.dumps(doc_list_one[ac][k], indent=4)
                        #print json.dumps(doc_list_two[ac][k], indent=4)
                        print ac, doc_list_one[ac][k], doc_list_two[ac][k]
                        continue

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

