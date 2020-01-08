graphIri="http://sparql.glygen.org"


#isql 1111 dba dba exec="SPARQL SELECT count(*) FROM <$graphIri> WHERE {?s ?p ?o .};"

#isql 1111 dba dba exec="SPARQL SELECT ?gsite_uri WHERE { ?gsite_uri <http://purl.jp/bio/12/glyco/conjugate#has_saccharide> ?glycan_uri . };"

isql 1111 dba dba exec="SPARQL SELECT ?s ?p ?o FROM <$graphIri> WHERE { ?s ?p ?o . };"
