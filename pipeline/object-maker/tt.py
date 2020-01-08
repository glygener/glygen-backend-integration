from Bio import AlignIO
from Bio.Align import AlignInfo

aln_file = "alignments/isoformset.uniprotkb.P14210-1.aln"
align = AlignIO.read(aln_file, "clustal")
#for record in align:
#    print record.id

summary_align = AlignInfo.SummaryInfo(align)
consensus = summary_align.dumb_consensus(1.0, ":")
print consensus

print summary_align

