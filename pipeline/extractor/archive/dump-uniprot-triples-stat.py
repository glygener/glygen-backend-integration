import os,sys
import json
import csv
from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--configfile",action="store",dest="configfile",help="config JSON file")
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

    (options,args) = parser.parse_args()
    for file in ([options.configfile, options.species]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    config_obj = json.loads(open(options.configfile, "r").read())
    species = options.species
    tax_id = "9606" if species == "human" else "10090"


    sparql = SPARQLWrapper("http://localhost:8890/sparql")

    sparql.setReturnFormat(JSON)
    graph_uri = "http://purl.glygen.org#uniprot_%s" % (species)

    query_header = """
        PREFIX faldo: <http://biohackathon.org/resource/faldo#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
    """
    
    class_query = query_header + " SELECT count(distinct ?s) as ?n FROM <%s> WHERE { ?s rdf:type %s . }"
    predicate_query = query_header + " SELECT count(distinct *) as ?n FROM <%s> WHERE { ?s %s ?o . }"
    prdval_query = query_header +" SELECT count(distinct *) as ?n FROM <%s> WHERE { ?s %s %s . }"
    clsprd_query = query_header + " SELECT count(distinct *) as ?n FROM <%s> WHERE { ?s rdf:type %s . ?s %s ?o . }"
    clsprdval_query = query_header + " SELECT count(distinct *) as ?n FROM <%s> WHERE { ?s rdf:type %s . ?s %s %s . }"


    out_json = {"class":{}, "predicate":{}, "predicatevalue":{},"classpredicate":{}, "classpredicatevalue":{}}

    for cls in config_obj["class"]:
        if cls.find("*") != -1:
            continue
        qs = class_query % (graph_uri, cls) 

        sparql.setQuery(qs)

        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            out_json["class"][cls] = int(result["n"]["value"])


    for prd in config_obj["predicate"]:
        if prd.find("*") != -1:
            continue
        qs = predicate_query % (graph_uri, prd)
        sparql.setQuery(qs)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            out_json["predicate"][prd] = int(result["n"]["value"])


      
    for prdval in config_obj["predicatevalue"]:
        prd_ns,prd_name,val_ns,val_type,val = prdval.split(":")
        prd = "%s:%s" %(prd_ns, prd_name)
        val_uri = "\"%s\"^^<%s%s>" % (val, config_obj["nsmap"][val_ns], val_type)
        qs = prdval_query % (graph_uri, prd, val_uri)
        sparql.setQuery(qs)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            out_json["predicatevalue"][prdval] = int(result["n"]["value"])

                
    for clsprd in config_obj["classpredicate"]:
        cls_ns,cls_name,prd_ns,prd_name = clsprd.split(":")
        cls = "%s:%s" %(cls_ns, cls_name)
        prd = "%s:%s" %(prd_ns, prd_name)
        qs = clsprd_query % (graph_uri,cls,prd)            
        sparql.setQuery(qs)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            out_json["classpredicate"][clsprd] = int(result["n"]["value"])


    for clsprdval in config_obj["classpredicatevalue"]:
        cls_ns,cls_name,prd_ns,prd_name,val_ns,val_type,val = clsprdval.split(":")
        cls = "%s:%s" %(cls_ns, cls_name)
        prd = "%s:%s" %(prd_ns, prd_name)
        val_uri = "\"%s\"^^<%s%s>" % (val, config_obj["nsmap"][val_ns], val_type)
        qs = clsprdval_query % (graph_uri,cls,prd,val_uri)
        sparql.setQuery(qs)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            out_json["classpredicatevalue"][clsprdval] = int(result["n"]["value"])

    print json.dumps(out_json, indent=4, sort_keys=True)

if __name__ == '__main__':
	main()
