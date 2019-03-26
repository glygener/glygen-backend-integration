
#dirName="/data/projects/glygen/generated/uniprot/n3/homo-sapiens/"
dirName="/data/projects/glygen/generated/uniprot/n3/mus-musculus/"

graphIri="http://glygendata.com#"

#isql 1111 dba dba exec="SPARQL CLEAR GRAPH <$graphIri>;"
isql 1111 dba dba exec="delete from DB.DBA.load_list;"
#isql 1111 dba dba exec="LOAD rdfloader.sql;"
isql 1111 dba dba exec="ld_dir ('$dirName', '*.n3', '$graphIri');"
isql 1111 dba dba exec="rdf_loader_run();"
isql 1111 dba dba exec="rdf_loader_run();"

wait

isql 1111 dba dba exec="checkpoint;"


#isql 1111 dba dba exec="SPARQL SELECT * FROM <$graphIri> WHERE {?s ?p ?o};"


