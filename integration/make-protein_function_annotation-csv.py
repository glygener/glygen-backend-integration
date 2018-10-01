import os,sys
import util
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


__version__="1.0"
__status__ = "Dev"


###############################
def get_refseq_dataframe(in_file, map_file):

    

    refseq2canonical = {}
    with open(map_file, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_frame:
            row_count += 1
            refseq2canonical[row[3]] = row[0]
        
    data_frame = {}
    field = "xxx"
    with open(in_file, "r") as FR:
        for line in FR:
            newfield = line[0:12].strip() 
            field = newfield if len(newfield) > 0 else field
            value = line[12:].strip()
            if field == "VERSION":
                ac = line.split(" ")[-1].strip()
                data_frame[ac] = {"ac":ac, "references":[]}
            elif field == "REFERENCE":
                data_frame[ac]["references"].append({"pubmed":"", "remark":""})
            elif field == "REMARK":
                data_frame[ac]["references"][-1]["remark"] += value.strip() + " "
            elif field == "PUBMED":
                data_frame[ac]["references"][-1]["pubmed"] = value

    rows = []
    for ac in data_frame:
        for o in data_frame[ac]["references"]:
            if o["remark"].find("GeneRIF:") != -1:
                o["remark"] = o["remark"][8:].strip()
                o["remark"] = o["remark"][0].upper() + o["remark"][1:]
                o["remark"] = o["remark"].replace("\"", "`")
                uniprot_canonical_ac = "N/A"
                if ac in refseq2canonical:
                    uniprot_canonical_ac = refseq2canonical[ac]
                    row = [uniprot_canonical_ac, "RefSeq",ac, o["pubmed"], o["remark"]]
                    rows.append(row)
    return rows




###############################
def get_uniprot_dataframe(in_file):



    data_frame = {}
    with open(in_file, 'r') as FR:
        csv_frame = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in csv_frame:
            row_count += 1
            if row_count == 1:
                continue
            data_frame[row[0]] = {"ac":row[0], "annotation":row[-1]}
       
    rows = []
    for canon in data_frame:
        if len(data_frame[canon]["annotation"]) > 0:
            word_list = data_frame[canon]["annotation"].split(" ")
            pmid_list = []
            for word in word_list:
                if word.find("PubMed:") != -1:
                    word = word.split(")")[0]
                    word = word.replace("(", "").replace(")", "").replace(",", "")
                    pmid_list.append(word[7:])
            evidence = "|".join(sorted(set(pmid_list)))
            row = [canon, "UniProtKB", canon, evidence, data_frame[canon]["annotation"]]
            rows.append(row)

    return rows



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--configfile",action="store",dest="configfile",help="config JSON file")
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

    (options,args) = parser.parse_args()
    for file in ([options.species, options.configfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    config_obj = json.loads(open(options.configfile, "r").read())
    species = options.species
    tax_id = "9606" if species == "human" else "10090"

    in_file_one = config_obj["pathinfo"]["downloads"] + "/refseq/protein/refseq_protein_all_%s.gp" % (tax_id)
    in_file_two = config_obj["pathinfo"]["reviewed"] + "/%s_protein_function.csv" % (species)
    map_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_refseq_mapping.csv" % (species)

    refseq_rows = get_refseq_dataframe(in_file_one, map_file)
    uniprot_rows = get_uniprot_dataframe(in_file_two)

    row = ["uniprotkb_acc_canonical", "database", "database_id", "evidence","annotation"]
    print "\"%s\"" % ("\",\"".join(row))
    for row in uniprot_rows + refseq_rows:
        print "\"%s\"" % ("\",\"".join(row))




if __name__ == '__main__':
	main()


