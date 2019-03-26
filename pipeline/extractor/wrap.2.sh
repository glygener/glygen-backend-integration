for ds in xref_reactome xref_hgnc xref_mgi xref_kegg xref_bgee xref_biomuta xref_brenda xref_cazy xref_cdd xref_chembl xref_dbsnp xref_ensembl xref_enzyme xref_genecards xref_geneid xref_interpro xref_mim xref_nextprot xref_oma xref_orthodb xref_panther xref_pdb xref_pfam xref_pro xref_refseq xref_unicarbkb; do
    python extract-from-endpoint.py  -s human -d $ds  > unreviewed/human_protein_"$ds".csv
    python extract-from-endpoint.py  -s mouse -d $ds  > unreviewed/mouse_protein_"$ds".csv
done
