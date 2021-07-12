import os
import csv
import sys
import json
import commands
import glob
import gzip
from optparse import OptionParser

sys.path.append('../../glytools/')
import libgly



def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version=" ")
    parser.add_option("-b","--batch",action="store",dest="batch",help="1, 2 ...")

    (options,args) = parser.parse_args()
    for file in ([options.batch]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    batch = int(options.batch)

    config_obj = json.loads(open("conf/config.json", "r").read())
 
    global path_obj
    path_obj = config_obj["pathinfo"]

    cid2glytoucan = {}
    data_frame = {}
    in_file = "unreviewed/glycan_xref_pubchem.csv"

    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "glycan_xref_pubchem_compound":
            cid2glytoucan[xref_id] = ac

    batch_size = 10

    start = (batch-1) * batch_size
    end = start + batch_size
    file_list = sorted(glob.glob("downloads/pubchem/compound/current/sdf4glygen/*.sdf.gz"))

    out_file = "generated/pubchem/compound/cid2inchi.batch.%s.csv" % (batch)
    FW = open(out_file, "w")
    FW.close()

    log_file = "logs/pubchem-downlads-step3.%s.log" % (batch)
    with open(log_file, "w") as FL:
        FL.write("started logging\n")


    seen = {}
    for in_file in file_list[start:end]:
        file_name = in_file.split("/")[-1]
        with open(log_file, "a") as FL:
            FL.write("started processing %s\n" % (file_name))
        prev_line = ""
        cid_count = 0
        with gzip.open(in_file, 'rb') as FR:
            cid = ""
            for line in FR:
                line = line.strip()
                if prev_line == "> <PUBCHEM_COMPOUND_CID>":
                    cid = line.strip()
                    newrow = [line]
                elif prev_line == "> <PUBCHEM_IUPAC_INCHI>":
                    newrow += [line]
                elif prev_line == "> <PUBCHEM_IUPAC_INCHIKEY>":
                    newrow += [line]
                    if cid not in seen and cid in cid2glytoucan and len(newrow) == 3:
                        with open(out_file, "a") as FA:
                            FA.write("\"%s\"\n" % ("\",\"".join(newrow)))
                        cid_count += 1
                        seen[cid] = True
                prev_line = line
        with open(log_file, "a") as FL: 
            FL.write("finished processing %s, mapped %s CIDs\n" % (file_name, cid_count))

                


if __name__ == '__main__':
        main()
