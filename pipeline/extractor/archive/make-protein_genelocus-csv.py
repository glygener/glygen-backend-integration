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


##########################
def load_isoform2locusinfo(ep_data):

    #?trsid gly:reverseStrand ?strand .
    #?trsid gly:transcriptRange ?rangeuri .

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    query_string = """
        PREFIX faldo: <http://biohackathon.org/resource/faldo#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?isoformuri ?trsid ?pepid ?chrid ?startpos ?endpos
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?isoformuri gly:ensTranscript ?trsid .
            ?trsid up:translatedTo ?pepid . 
            ?trsid gly:chromosome ?chrid .
            ?trsid gly:transcriptRange ?rangeuri .
            ?rangeuri faldo:begin ?beginuri .
            ?rangeuri faldo:end ?enduri .
            ?beginuri faldo:position ?startpos .
            ?enduri faldo:position ?endpos .
        }
        """


    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        isoform = result["isoformuri"]["value"].split("/")[-1]
        if isoform not in ep_data["isoform2locusinfo"]:
            ep_data["isoform2locusinfo"][isoform] = {
                    "chrid":result["chrid"]["value"].strip()
                    ,"trsid":result["trsid"]["value"].split("/")[-1] 
                    ,"pepid":result["pepid"]["value"].split("/")[-1]
                    ,"strand":"x"
                    ,"startpos":result["startpos"]["value"]
                    ,"endpos":result["endpos"]["value"]
            }

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


    ep_data = {
        "isoform2locusinfo":{}, 
    }
    load_isoform2locusinfo(ep_data)


    data_frame = {}
    in_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_idmapping.csv" % (species)
    load_dataframe(data_frame, "idmapping", in_file, ",")
  
    row = ["uniprot_canonical_ac","uniprot_isoform_ac","transcript_id","peptide_id","chr_id","start_pos","end_pos"]
    print "\"%s\"" % ("\",\"".join(row))
    for canon in data_frame["idmapping"]:
        for obj in data_frame["idmapping"][canon]:
            isoform_list = obj["reviewed_isoforms"] + obj["unreviewed_isoforms"]
            for isoform in isoform_list:
                if isoform in ep_data["isoform2locusinfo"]:
                    o = ep_data["isoform2locusinfo"][isoform]
                    row = [canon,isoform,o["trsid"], o["pepid"],
                            str(o["chrid"]),str(o["startpos"]),str(o["endpos"])]
                    print "\"%s\"" % ("\",\"".join(row))










if __name__ == '__main__':
	main()
