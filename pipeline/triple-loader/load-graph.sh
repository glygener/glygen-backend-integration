#species="human"
species="mouse"
#species="all"

#graphIri="http://sparql.glygen.org#uniprot_$species"
#dirName="/data/projects/glygen/generated/sparql/$species/"

#graphIri="http://sparql.glygen.org#disease"
#dirName="/data/projects/glygen/generated/sparql/disease/"

graphIri="http://sparql.glygen.org"
dirName="/data/projects/glygen/generated/sparql/all/"


isql 1111 dba dba exec="delete from DB.DBA.load_list;"
#isql 1111 dba dba exec="LOAD rdfloader.sql;"
isql 1111 dba dba exec="ld_dir ('$dirName', '*.nt', '$graphIri');"
isql 1111 dba dba exec="rdf_loader_run();"
isql 1111 dba dba exec="rdf_loader_run();"

wait

isql 1111 dba dba exec="checkpoint;"




