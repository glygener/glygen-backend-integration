graphIri="http://purl.glygen.org#uniprot_human"


isql 1111 dba dba exec="SPARQL SELECT count(distinct ?s) as ?n FROM <$graphIri> WHERE { ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://purl.uniprot.org/core/Protein/> . }"


