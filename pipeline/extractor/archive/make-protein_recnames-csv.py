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

#######################
def load_ac2fullname(ep_data, tax_id):

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
        
    sparql.setQuery("""
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX gly: <http://purl.glygen.org/>
        PREFIX up: <http://purl.uniprot.org/core/>
        SELECT ?ac ?fullname
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:recommendedName ?nameuri .
            ?nameuri up:fullName ?fullname . 
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
        }
    """)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        ep_data["ac2fullname"][ac] = result["fullname"]["value"].strip()
    return


#######################
def load_ac2shortname(ep_data, tax_id):

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
        
    sparql.setQuery("""
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX gly: <http://purl.glygen.org/>
        PREFIX up: <http://purl.uniprot.org/core/>
        SELECT ?ac ?shortname
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:recommendedName ?nameuri .
            ?nameuri up:shortName ?shortname .
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
        }
    """)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        ep_data["ac2shortname"][ac] = result["shortname"]["value"].strip()
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
        "ac2fullname":{}, 
        "ac2shortname":{},
    }
    load_ac2fullname(ep_data,tax_id)
    load_ac2shortname(ep_data,tax_id)

    data_frame = {}
    in_file = config_obj["pathinfo"]["reviewed"] + "/%s_protein_idmapping.csv" % (species)
    load_dataframe(data_frame, "idmapping", in_file, ",")
  
    row = ["uniprotkb_acc_canonical","recommended_name_full","recommended_name_short"]
    print "\"%s\"" % ("\",\"".join(row))
    for canon in data_frame["idmapping"]:
        ac = canon.split("-")[0]
        full_name = ep_data["ac2fullname"][ac] if ac in ep_data["ac2fullname"] else ""
        short_name = ep_data["ac2shortname"][ac] if ac in ep_data["ac2shortname"] else ""
        row = [canon, full_name, short_name]
        print "\"%s\"" % ("\",\"".join(row))
            








if __name__ == '__main__':
	main()
