import os
import csv
import sys
import json
import commands
import glob
import gzip


sys.path.append('../../glytools/')
import libgly



def remove_empty_files():

    cmd = "ls -l %spubchem/compound/sdf/*.sdf.gz" % (path_obj["downloads"])
    lines = commands.getoutput(cmd).split("\n")
    for line in lines:
        file_size = line.split()[4]
        out_file = line.split()[-1]
        if file_size == "0":
            cmd = "rm %s" % (out_file)
            x = commands.getoutput(cmd)




def main():


    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]
 
    global path_obj
    path_obj = config_obj["pathinfo"]


    cid2glytoucan = {}
    data_frame = {}
    for in_file in glob.glob("unreviewed/glycan_xref_pubchem.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[f_list.index("glytoucan_ac")]
            database_id = row[f_list.index("database_id")]
            if database_id[0:3] == "CID":
                cid = database_id[3:]
                cid2glytoucan[cid] = ac
                

    
    FL = open("logs/cid2inchi.log", "w")
    FL.close()

    FW = open("generated/pubchem/compound/cid2inchi.csv", "w")
    newrow = ["pubchem_cid", "inchi", "inchikey"]
    FW.write("\"%s\"\n" % ("\",\"".join(newrow)))

    remove_empty_files()

    seen = {}
    url_tmplt = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_%s_%s.sdf.gz"
    start = 1
    end = 25000
    n_failed = 0
    for i in xrange(0, 10000):
        start = i*25000 + 1
        end = start + 25000 - 1
        start = "000000000"[0:-len(str(start))] + str(start)
        end = "000000000"[0:-len(str(end))] + str(end)
        url = url_tmplt % (start, end)
        out_file = path_obj["downloads"] + "pubchem/compound/sdf/Compound_%s_%s.sdf.gz" % (start, end)
        if os.path.isfile(out_file) == False:
            cmd = "wget -O %s %s" % (out_file, url)
            x = commands.getoutput(cmd)
            remove_empty_files()

        if os.path.isfile(out_file) == True:
            prev_line = ""
            with gzip.open(out_file, 'rb') as FR:
                cid = ""
                for line in FR:
                    line = line.strip()
                    if prev_line == "> <PUBCHEM_COMPOUND_CID>":
                        cid = line
                        newrow = [line]
                    elif prev_line == "> <PUBCHEM_IUPAC_INCHI>":
                        newrow += [line]
                    elif prev_line == "> <PUBCHEM_IUPAC_INCHIKEY>":
                        newrow += [line]
                        if cid not in seen and cid in cid2glytoucan and len(newrow) == 3:
                            FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
                            seen[cid] = True

                    prev_line = line
        with open("logs/cid2inchi.log", "a") as FL:
            FL.write("Finished downloading and processing %s\n" % (out_file) )

    FW.close()

                


if __name__ == '__main__':
        main()
