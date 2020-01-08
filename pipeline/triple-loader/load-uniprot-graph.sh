# Make sure ../triple-maker/wrap-make-triples.sh is run
# to generate files under generated/sparql/glygen/


if [ -z "$1" ]
then
    echo ""
    echo "no species [human|mouse|rat|all] speficified!"
    echo ""
    exit
fi

species=$1

declare -A pattern_dict
pattern_dict["human"]="uniprot-proteome-homo-sapiens.*.nt"
pattern_dict["mouse"]="uniprot-proteome-mus-musculus.*.nt"
pattern_dict["rat"]="uniprot-proteome-rattus-norvegicus.*.nt"
pattern_dict["all"]="uniprot-proteome-*.nt"

glob_pattern=${pattern_dict[$species]}

graphIri="http://sparql.glygen.org#uniprot_$species"
dirName="/data/projects/glygen/generated/sparql/glygen/"


#Now load triples
isql 1111 dba dba exec="delete from DB.DBA.load_list;"
isql 1111 dba dba exec="ld_dir ('$dirName', '$glob_pattern', '$graphIri');"
isql 1111 dba dba exec="rdf_loader_run();"
isql 1111 dba dba exec="rdf_loader_run();"

wait

isql 1111 dba dba exec="checkpoint;"




