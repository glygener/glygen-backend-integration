species="human"
#species="mouse"
#species="all"

#graphIri="http://sparql.glygen.org#uniprot_$species"
graphIri="http://sparql.glygen.org"


isql 1111 dba dba exec="SPARQL SELECT count(*) FROM <$graphIri> WHERE {?s ?p ?o .};"

