graphIri="http://sparql.glygen.org#disease"

isql 1111 dba dba exec="SPARQL SELECT count(*) FROM <$graphIri> WHERE {?s ?p ?o .};"

