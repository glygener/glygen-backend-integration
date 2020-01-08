if [ -z "$1" ]
then
    echo ""
    echo "no species [human|mouse|rat|all] speficified!"
    echo ""
    exit
fi

species=$1
graphIri="http://sparql.glygen.org#uniprot_$species"


isql 1111 dba dba exec="SPARQL SELECT count(*) FROM <$graphIri> WHERE {?s ?p ?o .};"

