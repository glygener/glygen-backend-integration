ebi_dir="downloads/ebi/current/"
out_dir="generated/sparql/glygen/"

#Clear all nt files in out_dir
rm $out_dir/*.nt

#Make glygen triples
python make-glycan-triples.py  > $out_dir/glycan.nt
python make-proteoform-triples.py > $out_dir/proteoform.nt


#Partition EBI triples
python partition-nt.py -i $ebi_dir/uniprot-proteome-homo-sapiens.nt -o $out_dir
python partition-nt.py -i $ebi_dir/uniprot-proteome-mus-musculus.nt -o $out_dir
python partition-nt.py -i $ebi_dir/uniprot-proteome-rattus-norvegicus.nt -o $out_dir



