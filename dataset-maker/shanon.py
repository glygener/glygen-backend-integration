from Bio import AlignIO
from scipy.stats import entropy
import numpy


# given col number, and list of amino acids to consider
col = 1
aa_list = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "-"]

msa = AlignIO.read("tmp/example.aln", "clustal")

#get list of amino acids observed at given column
aa_list_at_col = []
for record in msa:
    aa = record.seq[col-1]
    aa_list_at_col.append(aa)

#dictionary of background counts
background_dict = {
    "A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,
    "I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,
    "R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"-":0
}

#calculate counts 
count_list = []
for aa in aa_list:
    freq_at_col = aa_list_at_col.count(aa)
    freq_background = background_dict[aa]
    count_list.append(freq_at_col + freq_background)

dist = list(count_list/numpy.linalg.norm(count_list))

#calculate entropy
shanon_ent = entropy(dist, base=2)
print shanon_ent
