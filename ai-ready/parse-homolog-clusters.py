import os
import sys
import json
import glob
from optparse import OptionParser


import libgly


def load_glycodict(gtype_list):

    file_list = glob.glob("unreviewed/*_proteoform_glycosylation_sites_*.csv")
    out_dict = {}
    gtc2species = {}
    seen = {}
    gtc2species = {}
    for in_file in file_list:
        species = in_file.split("/")[-1].split("_")[0]
        data_frame = {}
        delim = "," if in_file.split(".")[-1] == "csv" else "\t"
        libgly.load_sheet(data_frame, in_file, delim)
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            glytoucan_ac = row[f_list.index("saccharide")]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
            g_type = row[f_list.index("glycosylation_type")]
            if g_type not in gtype_list:
                continue
            if xref_key != "protein_xref_pubmed":
                continue
            if canon == "" or glytoucan_ac == "" or aa_pos == "":
                continue
            if canon not in out_dict:
                out_dict[canon] = {}
            combo_id = "%s|%s|%s" % (canon, aa_pos, glytoucan_ac)
            aa_pos = int(aa_pos)
            if aa_pos not in out_dict[canon]:
                out_dict[canon][aa_pos] = {}
            out_dict[canon][aa_pos][glytoucan_ac] = True
            
            if combo_id not in seen and glytoucan_ac != "":
                if glytoucan_ac not in gtc2species:
                    gtc2species[glytoucan_ac] = {}
                if species not in gtc2species[glytoucan_ac]:
                    gtc2species[glytoucan_ac][species] = 0
                gtc2species[glytoucan_ac][species] += 1
            seen[combo_id] = True

    return out_dict, gtc2species



def main():

    gtype_list = ["N-linked"]
    glyco_dict, gtc2species = load_glycodict(gtype_list)


    alphabet = "ACDEFGHIKLMNPQRSTVWYXBZJ"
    rel_dir = "/data/shared/glygen/releases/data/v-2.0.3/"
    file_list = glob.glob(rel_dir + "jsondb/alignmentdb/homologset.oma.*.json")
    #file_list = glob.glob(rel_dir + "jsondb/alignmentdb/homologset.mgi.41488513.json")

    for in_file in file_list:
        doc = json.loads(open(in_file, "r").read())
        cls_id = doc["cls_id"]
        alnpos2seqpos = {}
        canon2taxid = {}
        for obj in doc["sequences"]:
            canon, tax_id = obj["id"], obj["tax_id"]
            canon2taxid[canon] = tax_id
            seq_pos = 0
            for aln_pos in range(0, len(obj["aln"])):
                if alphabet.find(obj["aln"][aln_pos].upper()) != -1:
                    seq_pos += 1
                if aln_pos not in alnpos2seqpos:
                    alnpos2seqpos[aln_pos] = {}
                alnpos2seqpos[aln_pos][canon] = seq_pos
        
        for aln_pos in alnpos2seqpos:
            seen_taxid = {}
            canon2glycanlist  = {}
            for canon in alnpos2seqpos[aln_pos]:
                tax_id = canon2taxid[canon]
                seq_pos = alnpos2seqpos[aln_pos][canon] 
                glycan_list = []  
                if canon in glyco_dict:
                    if seq_pos in glyco_dict[canon]:
                        canon2glycanlist[canon] = {"tax_id":tax_id, "glycans":list(glyco_dict[canon][seq_pos].keys())}
                        seen_taxid[tax_id] = True
            cls_size = len(list(alnpos2seqpos[aln_pos].keys()))
            tax_idlist = list(set(seen_taxid.keys()))
            if len(tax_idlist) > 1:
                o = {"cls_id":cls_id, "cls_size":cls_size, "aln_pos":aln_pos, 
                        "glycosylation":canon2glycanlist}
                out_file = "temp/%s.%s.json" % (cls_id, aln_pos)
                with open(out_file, "w") as FW:
                    FW.write("%s\n" % (json.dumps(o, indent=4)))
            
    return


                


if __name__ == '__main__':
        main()
