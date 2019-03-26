import os,sys
import json
import csv

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 
import time
import datetime



__version__="1.0"
__status__ = "Dev"


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

    species_short = options.species
    species_dict = {
            "human":"homo sapiens"
            ,"mouse":"mus musculus" 
    }
    species_name = species_dict[species_short]


    sparql = SPARQLWrapper("http://203.101.226.16:40935/unicarbkbv2.0.1/query")
    query = '''
	prefix gc: <http://www.oegov.org/core/owl/gc#>
	prefix sio: <http://semanticscience.org/resource/>
	prefix gl: <http://schema.geolink.org/>
	prefix owl: <http://www.w3.org/2002/07/owl#> 
	prefix ms: <http://www.berkeleybop.org/ontologies/ms.owl> 
	prefix xsd: <http://www.w3.org/2001/XMLSchema#> 
	prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
	prefix gco: <http://purl.jp/bio/12/glyco/conjugate#> 
	prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
	prefix unicarbkb: <http://purl.org/unicarbkb> 
	prefix faldo: <http://www.biohackathon.org/resource/faldo/> 
	prefix pato: <http://purl.obolibrary.org/obo/uo.owl>
	prefix dcterms: <http://purl.org/dc/terms/> 
	prefix bibo: <http://purl.org/ontology/bibo/> 
	prefix unicorn: <http://purl.jp/bio/12/glyco/unicorn/> 
	prefix uniprot: <http://purl.uniprot.org/core/> 
	prefix foaf: <http://xmlns.com/foaf/0.1/> 
	prefix dc: <http://purl.org/dc/elements/1.1/> 
	prefix glycan: <http://purl.jp/bio/12/glyco/glycan/>
	
        SELECT distinct ?Protein ?Position ?Pmid ?Id ?TypeAminoAcid
        WHERE {
          ?ReferencedCompound gco:ReferencedProtein ?Protein ;
           glycan:published_in ?Published ; glycan:is_from_source ?Source .
          ?Published glycan:has_pmid ?Pmid .
 	  ?Protein gco:has_protein ?ProteinAcc .
          ?Source glycan:has_taxon ?Taxon .
          ?Taxon uniprot:scientificName ?Species .
          ?Protein gco:glycosylated_at ?Region .
          ?Region faldo:ExactPosition ?Faldo .
          ?Faldo faldo:position ?Position .
          ?Region gco:has_saccharide_set ?Set .
          ?Set sio:is-component-part-of ?SetItem .
          ?SetItem owl:sameAs ?Saccharide .
          ?Saccharide dcterms:identifier ?Id .
          ?Faldo gco:has_amino_acid ?AminoAcid .
          ?AminoAcid gco:amino_acid ?TypeAminoAcid .

          FILTER regex(str(?Species), "%s") .
        }'''

    query = query % (species_name)

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
  
    time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')
    out_file = config_obj["pathinfo"]["downloads"] + "/uckb/%s_glycosylation_%s.csv" % (species_short,time_stamp)
    FW = open(out_file, "w")
    for result in results["results"]["bindings"]:
        FW.write("%s,%s,%s,%s,%s\n" % (
            result["Protein"]["value"], 
            result["Position"]["value"],
            result["Pmid"]["value"],
            result["Id"]["value"],
            result["TypeAminoAcid"]["value"]
        ))
    FW.close()




if __name__ == '__main__':
	main()
