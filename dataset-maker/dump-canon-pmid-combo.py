import os,sys
import libgly3
import glob
import subprocess
from optparse import OptionParser




def load_ac2canon(species):

    dict_one, dict_two = {}, {}
    data_frame = {}
    in_file = reviewed_dir + "%s_protein_masterlist.csv" % (species)
    libgly3.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        dict_one[canon] = canon
        isoform_list = [row[f_list.index("reviewed_isoforms")],
                        row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            dict_one[ac] = canon
            dict_one[isoform] = canon
            if isoform == canon:
                dict_two[ac] = canon

    return dict_one, dict_two



def main():


    
    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog ")
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

    (options,args) = parser.parse_args()
    for file in ([options.species]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    global reviewed_dir

    species = options.species
    download_dir = "2024_02_13"
    download_dir = "2024_04_08"
    reviewed_dir = "unreviewed/"
    #reviewed_dir = "/data/shared/glygen/releases/data/v-2.4.1/reviewed/"



    nt_file = ""
    cmd = "grep %s generated/misc/species_info.csv" % (species)
    line_list = subprocess.getoutput(cmd).split("\n")
    for line in line_list:
        row = line.split(",")
        if row[-3].find(".nt") != -1:
            nt_file = "downloads/ebi/%s/%s" % (download_dir, row[-3])

    
    seen = {}
    ac2canon, ac2canon_strict = load_ac2canon(species)
    cmd = "grep \" <http://purl.uniprot.org/core/citation> \" " + nt_file
    line_list = subprocess.getoutput(cmd).split("\n")
    for line in line_list:
        row = line.strip().split("/")
        ac = row[4].split(">")[0]
        pmid = row[12].split(">")[0] 
        if ac in ac2canon:
            canon = ac2canon[ac]
            combo = "%s|%s" % (canon, pmid)
            seen[combo] = True


    file_list = [reviewed_dir +  "%s_protein_function_uniprotkb.csv" % (species)]
    file_list += [reviewed_dir +  "%s_protein_ptm_annotation_uniprotkb.csv" % (species)]
    file_list += [reviewed_dir +  "%s_protein_site_annotation_uniprotkb.csv" % (species)]
    for in_file in file_list:
        data_frame = {}
        libgly3.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            if xref_key.find("xref_pubmed") != -1:
                combo = "%s|%s" % (canon, xref_id)
                seen[combo] = True

    n_one = len(seen.keys())


    seen = {}
    in_file = reviewed_dir + "%s_protein_citations_uniprotkb.csv" % (species) 
    data_frame = {}
    libgly3.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        if len(f_list) != len(row):
            print (f_list)
            print (row)
            print ()
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key.find("xref_pubmed") != -1:
            combo = "%s|%s" % (canon, xref_id)
            seen[combo] = True

    n_two = len(seen.keys())
    print ("expected around: %s, found=%s" % (n_one, n_two))


    return



if __name__ == '__main__':
    main()



