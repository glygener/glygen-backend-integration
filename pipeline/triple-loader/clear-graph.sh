#species="human"
#species="mouse"
species="all"

#graphIri="http://sparql.glygen.org#uniprot_$species"

graphIri="http://sparql.glygen.org"


#s="http://sparql.uniprot.org/uniprot/P27449"
#p="http://sparql.uniprot.org/core/organism"
#o="http://sparql.uniprot.org/taxonomy/9606"

#isql 1111 dba dba exec="SPARQL DELETE DATA FROM <$graphIri> {?ac <$p> <$o> .};"

isql 1111 dba dba exec="SPARQL CLEAR GRAPH <$graphIri>;" 
