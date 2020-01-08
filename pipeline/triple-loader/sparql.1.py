import os,sys
import json
import ssl

from SPARQLWrapper import SPARQLWrapper, JSON 

###############################
def main():

    if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
        ssl._create_default_https_context = ssl._create_unverified_context

    #sparql = SPARQLWrapper("http://localhost:8880/sparql")
    sparql = SPARQLWrapper("https://glygen.org:8890/sparql")
    #sparql = SPARQLWrapper("https://sparql.glygen.org:8890/sparql")

    sparql.setQuery("""
        SELECT ?s ?p ?o 
        FROM <http://sparql.glygen.org>
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
