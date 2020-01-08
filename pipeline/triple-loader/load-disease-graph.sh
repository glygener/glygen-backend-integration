graphIri="http://sparql.glygen.org#disease"
dirName="/data/projects/glygen/generated/sparql/disease/"

infile="/data/projects/glygen/downloads/do/current.owl"
outfile="$outDir/disease/disease.nt"

/usr/local/bin/rapper -o ntriples  $infile >  $outfile


#First clear graph
isql 1111 dba dba exec="SPARQL CLEAR GRAPH <$graphIri>;"

#Now load triples
isql 1111 dba dba exec="delete from DB.DBA.load_list;"
isql 1111 dba dba exec="ld_dir ('$dirName', '*.nt', '$graphIri');"
isql 1111 dba dba exec="rdf_loader_run();"
isql 1111 dba dba exec="rdf_loader_run();"

wait

isql 1111 dba dba exec="checkpoint;"




