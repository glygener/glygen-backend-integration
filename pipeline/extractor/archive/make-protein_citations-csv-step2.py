import os,sys
import json
import csv


from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


__version__="1.0"
__status__ = "Dev"


######################
def load_dataframe(data_frame, sheet_name, in_file, separator):

    data_frame[sheet_name] = {}
    with open(in_file, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=separator, quotechar='"')
        row_count = 0
        for row in csv_grid:
            row_count += 1
            if row_count == 1:
                field_list = row
            else:
                row_obj = {}
                for j in xrange(1,len(field_list)):
                    field_name  = field_list[j]
                    row_obj[field_name] = [] if row[j].strip() == "" else row[j].replace("\"","").split("|")
                main_id = row[0].strip()
                if main_id not in data_frame[sheet_name]:
                    data_frame[sheet_name][main_id] = []
                data_frame[sheet_name][main_id].append(row_obj)
    return




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


    data_frame = {}
    in_file = config_obj["pathinfo"]["intermediate"] + "/%s_protein_citations.csv" % (species)
    load_dataframe(data_frame, "citations", in_file, ",")

    in_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_blacklisted_pmids.csv" % (species)
    load_dataframe(data_frame, "pmidblacklist", in_file, ",")


    row = ["uniprot_canonical_ac","pmid","title","journal_name"]
    print "\"%s\"" % ("\",\"".join(row))
    for canon in data_frame["citations"]:
        rows = data_frame["citations"][canon]
        for row in rows:
            pmid_list = []
            for pmid in row["pmid"]:
                if pmid not in data_frame["pmidblacklist"]:
                    pmid_list.append(pmid)
            if len(pmid_list) > 0:
                row["pmid"] = "|".join(pmid_list)
                out_row = [canon,row["pmid"], row["title"][0], row["journal_name"][0]]
                print "\"%s\"" % ("\",\"".join(out_row))









if __name__ == '__main__':
	main()
