import os
import csv
import sys
import json
import commands
import glob
import gzip


sys.path.append('../../glytools/')
import libgly





def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]
 
    cid2glytoucan = {}
    data_frame = {}
    for in_file in glob.glob("unreviewed/*_glycan_xref_pubchem.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[f_list.index("glytoucan_ac")]
            database_id = row[f_list.index("database_id")]
            if database_id[0:3] == "CID":
                cid = database_id[3:]
                cid2glytoucan[cid] = ac


    seen = {}
    start = 1
    end = 25000
    n_failed = 0
    for i in xrange(0, 10000):
        start = i*25000 + 1
        end = start + 25000 - 1
        start = "000000000"[0:-len(str(start))] + str(start)
        end = "000000000"[0:-len(str(end))] + str(end)
        in_file = "downloads/pubchem/compound/sdf/Compound_%s_%s.sdf.gz" % (start, end)
        out_file = "downloads/pubchem/compound/sdf4glygen/Compound_%s_%s.sdf" % (start, end)
        if os.path.isfile(in_file) == True:
            FW = open(out_file, "w")
            flag = False 
            with gzip.open(in_file, 'rb') as FR:
                prev_line, cid, buf = "", "", ""
                for line in FR:
                    buf += line
                    line = line.strip()
                    if line == "$$$$":
                        if cid != "" and cid in cid2glytoucan:
                            flag = True
                            FW.write("%s" % (buf))
                        cid, buf = "", ""
                    if prev_line == "> <PUBCHEM_COMPOUND_CID>":
                        cid = line
                    prev_line = line
            FW.close()
            if flag == False:
                cmd = "rm -f %s" % (out_file)
                x = commands.getoutput(cmd)



if __name__ == '__main__':
        main()
