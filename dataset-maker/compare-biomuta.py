import os,sys
import json
import csv

from optparse import OptionParser

import commands
import glob
import csvutil


__version__="1.0"
__status__ = "Dev"



###############################
def main():


    map_file = "reviewed/human_protein_transcriptlocus.csv"
    old_in_file = "reviewed/human_protein_mutation_cancer.csv"
    #old_in_file = "downloads/biomuta/v-3.0/biomuta.csv"
    new_in_file = "downloads/biomuta/v-5.0/biomuta_tcga.csv"

    data_frame = {}
    csvutil.load_sheet(data_frame, map_file, [], ",")
    f_list = data_frame["fields"]
    pepid2ac = {}
    for row in data_frame["data"]:
        ac = row[0].split("-")[0]
        pep_id = row[3].strip()
        pepid2ac[pep_id] = ac


    seen = {"old":{}, "newone":{}, "newtwo":{}}
    in_new, mapped, unmapped, in_old, in_both = 0, 0, 0, 0, 0
    badrows = 0

    data_frame = {}
    csvutil.load_sheet(data_frame, new_in_file, [], ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        pep_id = row[f_list.index("ENSP")].strip()
        pos = row[f_list.index("Protein_position")]
        aa = row[f_list.index("Amino_acids")]
        if pos.find("/") == -1 or aa.find("/") == -1:
            badrows += 1
        else:
            aa_pos = pos.split("/")[0]
            ref_aa = aa.split("/")[0]
            alt_aa = aa.split("/")[1]
            tmp_list = [pep_id, aa_pos, ref_aa, alt_aa]
            mut = "|".join(tmp_list)
            if mut not in seen["newone"]:
                seen["newone"][mut] = True
                in_new += 1
                if pep_id not in pepid2ac:
                    unmapped += 1
                else:
                    mapped += 1
                    tmp_list = [pepid2ac[pep_id], aa_pos, ref_aa, alt_aa]
                    m = "|".join(tmp_list)
                    if m not in  seen["newtwo"]:
                        seen["newtwo"][m] = True


     
    data_frame = {}
    csvutil.load_sheet(data_frame, old_in_file, [], ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        source = row[f_list.index("source")] if "source" in f_list else ""
        source = row[f_list.index("data_source")] if "data_source" in f_list else source
        if source not in ["tcga"]:
            continue
        ac = ""
        if "uniprot_canonical_ac" in f_list:
            ac = row[f_list.index("uniprot_canonical_ac")].split("-")[0]
        if "uniprotkb_canonical_ac" in f_list:
            ac = row[f_list.index("uniprotkb_canonical_ac")].split("-")[0]

        tmp_list = [ac]
        for f in ["aa_pos","ref_aa","alt_aa"]:
            tmp_list.append(row[f_list.index(f)])
        mut = "|".join(tmp_list)
        if mut not in seen["old"]:
            in_old += 1
            seen["old"][mut] = True
            if mut in seen["newtwo"]:
                in_both += 1

    print "in_new=%s, mapped=%s, unmapped=%s" % (in_new,mapped,unmapped)
    print "in_old=%s" % (in_old)
    print "in_both=%s" % (in_both)
    print "badrows=%s" % (badrows)












if __name__ == '__main__':
        main()


