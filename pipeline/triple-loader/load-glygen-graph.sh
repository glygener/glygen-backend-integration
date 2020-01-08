# Make sure ../triple-maker/wrap-make-triples.sh is run
# to generate files under generated/sparql/glygen/


graphIri="http://sparql.glygen.org"
dirName="/data/projects/glygen/generated/sparql/glygen/"

isql 1111 dba dba exec="delete from DB.DBA.load_list;"
#isql 1111 dba dba exec="LOAD rdfloader.sql;"
isql 1111 dba dba exec="ld_dir ('$dirName', '*.nt', '$graphIri');"
isql 1111 dba dba exec="rdf_loader_run();"
isql 1111 dba dba exec="rdf_loader_run();"

wait

isql 1111 dba dba exec="checkpoint;"




