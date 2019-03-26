import os,sys
import json

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

import pymongo
from pymongo import MongoClient


__version__="1.0"
__status__ = "Dev"


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


    ac2canon = {}

    query_string = """
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
        SELECT ?ac ?isoformuri
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?ac up:sequence ?isoformuri .
            ?isoformuri gly:canonical 'true' ^^<http://www.w3.org/2001/XMLSchema#boolean> .
            ?ac up:organism <http://purl.uniprot.org/taxonomy/%s> .
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
        ac2canon[ac] = canon



    sparql.setQuery("""
        PREFIX rdf:   <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
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

    canon2fullname = {}
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        ac = result["ac"]["value"].split("/")[-1]
        canon = ac2canon[ac] if ac in ac2canon else ""
        if canon != "":
            full_name = result["fullname"]["value"].strip()
            canon2fullname[canon] = full_name


    client = MongoClient('mongodb://localhost:27017')
    db = client["glygen"]
    query_obj = {"species.taxid": {'$eq': int(tax_id)}}

    for doc in db["c_protein"].find(query_obj):
        canon = doc["uniprot_canonical_ac"]
        #if canon not in canon2fullname:
        #    print "A", canon
        if canon in canon2fullname:
            full_name = doc["recommendedname"]["full"]
            if full_name != canon2fullname[canon]:
                print "%s|%s|%s" % (canon, full_name, canon2fullname[canon])



if __name__ == '__main__':
	main()
