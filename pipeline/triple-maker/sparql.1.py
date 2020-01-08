import os,sys
import json

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


__version__="1.0"
__status__ = "Dev"


###############################
def main():

    
    #sparql = SPARQLWrapper("http://localhost:8890/sparql")
    #sparql = SPARQLWrapper("https://glygen.org:8890/sparql")
    sparql = SPARQLWrapper("https://sparql.glygen.org:8890/sparql")



    sparql.setQuery("""
        PREFIX gly: <http://purl.glygen.org/>
        SELECT ?s ?p ?o 
        FROM <http://purl.glygen.org#uniprot>
        WHERE { 
            ?s ?p ?o .
        }
        LIMIT 10
    """)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
        print result["s"]["value"], result["o"]["value"]



if __name__ == '__main__':
	main()
