import os,sys
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


######################
def load_dataframe(in_file, separator):

    data_frame = {}
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
                if main_id not in data_frame:
                    data_frame[main_id] = []
                data_frame[main_id].append(row_obj)

    return data_frame




###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

    (options,args) = parser.parse_args()
    for file in ([options.species]):
        if not (file):
            parser.print_help()
            sys.exit(0)
        
    species = options.species
    tax_id = "9606" if species == "human" else "10090"

    sparql = SPARQLWrapper("http://localhost:8890/sparql")


    seen = {}
    ac2canon = {}
    query_string = """
        PREFIX rdf:   <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?ac ?isoformuri
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:sequence ?isoformuri .
            ?isoformuri gly:canonical 'true' ^^<http://www.w3.org/2001/XMLSchema#boolean> .
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
        }
        """
    query_string = query_string % (tax_id)
    #print query_string
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        canon = result["isoformuri"]["value"].split("/")[-1]
        seen[canon] = True
        ac2canon[ac] = canon


    cls_dict = {"one":{}, "two":{}}

    query_string = """
        PREFIX rdf:   <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?ac ?isoformuri
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:sequence ?isoformuri .
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
            ?ac rdf:type <http://purl.uniprot.org/core/Protein> .
        }
        """
    query_string = query_string % (tax_id)
    #print query_string
    sparql.setQuery(query_string)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        isoform = result["isoformuri"]["value"].split("/")[-1]
        #canon = ac2canon[ac] if ac in ac2canon else "X_" + ac
        if ac in ac2canon:
            canon = ac2canon[ac]
            if canon not in cls_dict["two"]:
                cls_dict["two"][canon] = []
            cls_dict["two"][canon].append(isoform)

    n = 0
    for canon in cls_dict:
        if canon[0:2] != "X_":
            #print canon, ",".join(sorted(set(cls_dict[canon])))
            #print cls_dict[canon]
            n += len(cls_dict[canon])


    in_file = "reviewed/%s_protein_idmapping.csv" % (species)
    data_frame = load_dataframe(in_file, ",")
    for canon in data_frame:
        cls_dict["one"][canon] = []
        for obj in data_frame[canon]:
            cls_dict["one"][canon] += obj["reviewed_isoforms"] + obj["unreviewed_isoforms"]
    
        #if canon not in cls_dict["two"]:
        #    print "%s,%s" %(canon, "|".join(sorted(set(cls_dict["one"][canon]))))
        #if canon in cls_dict["one"] and canon in cls_dict["two"]:
        #    if sorted(set(cls_dict["one"][canon])) != sorted(set(cls_dict["two"][canon])):
        #        print "%s,%s,%s" %(canon, "|".join(sorted(set(cls_dict["one"][canon]))),"|".join(sorted(set(cls_dict["two"][canon]))))

    for canon in cls_dict["two"]:
        if canon not in cls_dict["one"]:
            print "%s,%s" %(canon, "|".join(sorted(set(cls_dict["two"][canon]))))




if __name__ == '__main__':
	main()
