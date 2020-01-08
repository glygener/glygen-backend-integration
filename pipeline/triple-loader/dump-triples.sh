graphIri="http://sparql.glygen.org"


#isql 1111 dba dba exec="SPARQL SELECT count(*) FROM <$graphIri> WHERE {?s ?p ?o .};"
isql 1111 dba dba exec="dump_one_graph ('http://sparql.glygen.org', './test_', 1000000000);"

