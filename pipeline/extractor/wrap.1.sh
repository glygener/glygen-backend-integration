
for ds in idmapping motif; do
    python expand-piped-values-type-two.py -i reviewed/human_glycan_"$ds".csv > unreviewed/human_glycan_"$ds".csv
done

for ds in idmapping motif; do
    python expand-piped-values-type-two.py -i reviewed/mouse_glycan_"$ds".csv > unreviewed/mouse_glycan_"$ds".csv
done


for ds in ovca_log all_glycosylation_sites glycosylation_statistics sequence_alignment; do
    python expand-piped-values-type-two.py -i reviewed/human_proteoform_"$ds".csv > unreviewed/human_proteoform_"$ds".csv
done

for ds in all_glycosylation_sites glycosylation_statistics; do
    python expand-piped-values-type-two.py -i reviewed/mouse_proteoform_"$ds".csv > unreviewed/mouse_proteoform_"$ds".csv

done


for ds in glycosylation_sites; do
    python expand-piped-values-type-two.py -i reviewed/hcv_proteoform_"$ds".csv > unreviewed/hcv_proteoform_"$ds".csv
done


