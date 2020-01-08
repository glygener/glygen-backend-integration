import os,sys

import json
import csv

from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


import commands
import glob
import sparqlutil


sys.path.append('../../glytools/')
import libgly


__version__="1.0"
__status__ = "Dev"



def get_blacklisted_pmids(species):

    black_list = []
    in_file = "compiled/%s_protein_blacklisted_pmids.csv" % (species)
    if os.path.isfile(in_file) == True:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        for row in data_frame["data"]:
            black_list.append(row[0])
        black_list = sorted(set(black_list))

    return black_list






def get_citation(pmid):

    row = []
    out_file = path_obj["downloads"] + "ncbi/pubmed/medline/pmid.%s.txt" % (pmid)
    if os.path.isfile(out_file) == True:
        obj = {}
        with open(out_file, "r") as FR:
            lcount = 0
            prev_key = ""
            for line in FR:
                lcount += 1
                if lcount > 3:
                    key = line[0:4].strip()
                    val = line[5:].strip()
                    if key not in obj:
                        obj[key] = []
                    if key == "":
                        obj[prev_key].append(val)
                    else:
                        obj[key].append(val)
                        prev_key = key
        if "TI" in obj and "JT" in obj and "DP" in obj and "AU" in obj:
            title = " ".join(obj["TI"]).replace("\"", "`")
            journal = " ".join(obj["JT"])
            pubdate = " ".join(obj["DP"])
            authors = ", ".join(obj["AU"])
            row = [pmid, title, journal, pubdate, authors]

    return row



def extract_citations_unicarbkb_ds(species):

    black_list = get_blacklisted_pmids(species)

    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_proteoform_glycosylation_sites_unicarbkb.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    newrow = ["uniprotkb_canonical_ac","pmid","title","journal_name","publication_date", "authors"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        pmid = row[f_list.index("evidence")]
        if pmid in black_list:
            continue
        combo_id = "%s %s" % (canon, pmid)
        newrow = get_citation(pmid)
        if newrow != []:
            if combo_id not in seen:
                print "\"%s\"" % ("\",\"".join([canon] + newrow))
            seen[combo_id] = True

    return



def extract_glycosylation_sites_glyconnect_ds(species):

    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon


    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_one2three = {}
    aa_three2one = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        aa_one2three[one] = three
        aa_three2one[three] = one

    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_allsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())


    in_file = path_obj["downloads"] + "glyconnect/current.json"
    obj_list = json.loads(open(in_file, "r").read())["glycosylations"]


    newrow = ["uniprotkb_canonical_ac","glycosylation_type","saccharide","amino_acid","glycosylation_site_uniprotkb","evidence", "composition","glyconnect_id","is_composition_only","glycan_core","glyconnect_protein_id","glytoucan_ac_composition"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for obj in obj_list:
        if obj["uniprot_acc"] not in ac2canon:
            continue
        canon = ac2canon[obj["uniprot_acc"]]
        if obj["taxonomy_id"] != str(species_obj[species]["taxid"]):
            continue
        
        aa_three = obj["glyco_site"].split("-")[0]
        if aa_three not in aa_three2one:
            continue
        

        aa_one = aa_three2one[aa_three]
        aa_pos = obj["glyco_site_location"]
        if aa_pos == 0 or aa_pos >= len(seq_hash[canon]):
            continue
        if aa_one != seq_hash[canon][aa_pos-1]:
            continue
        evidence = obj["pmid"] if "pmid" in obj else ""
        glyconnect_id = str(obj["structure_id"]) if "structure_id" in obj else ""
        glytoucan_ac = obj["structure_glytoucan_id"] if "structure_glytoucan_id" in obj else ""
        is_composition_only = obj["is_composition_only"] if "is_composition_only" in obj else ""
        glycan_core = obj["glycan_core"] if "glycan_core" in obj else ""
        glyconnect_protein_id = obj["protein_id"] if "protein_id" in obj else ""
        glytoucan_ac_composition = obj["composition_glytoucan_id"] if "composition_glytoucan_id" in obj else ""
        newrow = [canon,obj["glycan_type"],glytoucan_ac,aa_three,str(aa_pos),evidence,obj["composition"],glyconnect_id,is_composition_only,glycan_core,str(glyconnect_protein_id),glytoucan_ac_composition]
        print "\"%s\"" % ("\",\"".join(newrow))


    return


def extract_glycosylation_sites_literature_ds(species):

    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        ac2canon[canon] = canon

        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon
            ac2canon[isoform] = canon


    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_allsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())

    data_frame = {}
    in_file = path_obj["misc"] + "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_three2one = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        aa_three2one[three] = one


    newrow = ["uniprotkb_canonical_ac","glycosylation_type","saccharide","amino_acid",
            "glycosylation_site_uniprotkb","evidence"]
    print "\"%s\"" % ("\",\"".join(newrow))

    seen_row = {}
    log_file_one = path_obj["logs"] + "%s_proteoform_glycosylation_sites_literature.1.log" % (species)
    FL1 = open(log_file_one, "w")
    for in_file in glob.glob("compiled/glycosylation_sites_*"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            flag_dict = {
                "no_canon":False, "bad_pos":False, "aa_mismatch":False, "invalid_aa":False,
                "glycan_without_glytype":False,"bad_glytype":False,"no_pmid":False,
                "no_glycan_invalid_aa":False, "glycan_without_n_or_o_glytype": False
            }
            ac = row[f_list.index("uniprotkb_ac")]
            aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
            if row[f_list.index("glycosylation_site_uniprotkb")].isdigit() == False:
                flag_dict["bad_pos"] = True
            else:
                aa_pos = int(aa_pos)
            evidence = row[f_list.index("evidence")]
            aa_three = row[f_list.index("amino_acid")]
            glycosylation_type = row[f_list.index("glycosylation_type")]
            
            if ac not in ac2canon:
                flag_dict["no_canon"] = True
            canon = ac2canon[ac] if ac in ac2canon else ac
            seq = seq_hash[canon] if ac in ac2canon else ""

            if aa_three not in aa_three2one:
                flag_dict["invalid_aa"] = True
            aa_one = aa_three2one[aa_three] if aa_three in aa_three2one else ""
            if seq != "":
                if aa_pos == 0 or aa_pos >= len(seq):
                    flag_dict["bad_pos"] = True
                elif aa_one != seq[aa_pos-1]:
                    flag_dict["aa_mismatch"] = True
        
            newrow = [canon, glycosylation_type, "",aa_three, str(aa_pos), evidence]
            
            #To avoid redundancy
            row_str = "\",\"".join(newrow)
            if row_str in seen_row:
                continue
            seen_row[row_str] = True 
            
            flag_list = []
            for flag in flag_dict:
                if flag_dict[flag] == True:
                    flag_list.append(flag)
            if flag_list != []:
                newrow.append(";".join(flag_list))
                FL1.write("\"%s\"\n" % ("\",\"".join(newrow)))
            else:
                print "\"%s\"" % ("\",\"".join(newrow))

    FL1.close()


    return



def extract_glycosylation_sites_tyr_o_linked_ds(species):

    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon


    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_one2three = {}
    aa_three2one = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        aa_one2three[one] = three
        aa_three2one[three] = one

        
    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_allsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())

    data_frame = {}
    in_file = "compiled/%s_proteoform_glycosylation_sites_tyr_o_linked.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    newrow = ["uniprotkb_canonical_ac","glycosylation_type","saccharide","amino_acid","glycosylation_site_uniprotkb","evidence"]
    print "\"%s\"" % ("\",\"".join(newrow))

    log_file_one = path_obj["logs"] + "%s_proteoform_glycosylation_sites_tyr_o_linked.1.log" % (species)
    FL1 = open(log_file_one, "w")        
    seen = {}
    seen_row = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("uniprotkb_ac")]
        aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
        aa_three = row[f_list.index("amino_acid_uniprotkb")]
        glycosylation_type = row[f_list.index("glycosylation_type")].lower()
        evidence = row[f_list.index("evidence")]
        pmid = row[f_list.index("pmid")]
        flag_dict = {
            "no_canon":False, "bad_pos":False, "aa_mismatch":False, "invalid_aa":False,
            "glycan_without_glytype":False,"bad_glytype":False,"no_pmid":False,
            "no_glycan_invalid_aa":False, "glycan_without_n_or_o_glytype": False
        }
        if aa_pos.isdigit() == False:
            flag_dict["bad_pos"] = True
        else:
            aa_pos = int(aa_pos)
        if ac not in ac2canon:
            flag_dict["no_canon"] = True
        canon = ac2canon[ac] if ac in ac2canon else ac
        seq = seq_hash[canon] if ac in ac2canon else ""

        if aa_three not in aa_three2one:
            flag_dict["invalid_aa"] = True
        aa_one = aa_three2one[aa_three] if aa_three in aa_three2one else ""
        if seq != "":
            if aa_pos == 0 or aa_pos >= len(seq):
                flag_dict["bad_pos"] = True
            elif aa_one != seq[aa_pos-1]:
                flag_dict["aa_mismatch"] = True
        flag_list = []
        for flag in flag_dict:
            if flag_dict[flag] == True:
                flag_list.append(flag)

        newrow_list = [
            [canon, glycosylation_type, "",aa_three, str(aa_pos), evidence],
            [canon, glycosylation_type, "",aa_three, str(aa_pos), pmid]
        ]
        if flag_list != []:
            for newrow in newrow_list:
                newrowstr = json.dumps(newrow)
                if newrowstr not in seen_row:
                    newrow.append(";".join(flag_list))
                    FL1.write("\"%s\"\n" % ("\",\"".join(newrow)))
                    seen_row[newrowstr] = True
        else:
            for newrow in newrow_list:
                newrowstr = json.dumps(newrow)
                if newrowstr not in seen_row:
                    print "\"%s\"" % ("\",\"".join(newrow))
                    seen_row[newrowstr] = True



    FL1.close()

    return





def extract_glycosylation_sites_harvard_ds(species):


    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon


    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_one2three = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        aa_one2three[one] = three

    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_allsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())


    data_frame = {}
    in_file = path_obj["downloads"] + "harvard/%s_glycosylation_current.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    
    newrow = ["uniprotkb_canonical_ac","glycosylation_type","saccharide","amino_acid","glycosylation_site_uniprotkb","evidence", "composition"]
    print "\"%s\"" % ("\",\"".join(newrow))


    seen = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("UniProt Accession")]
        s_list = row[f_list.index("GlycoSite(s) (protein)")].split(";")
        g_list = row[f_list.index("Glycan composition")].split(";")
        evidence = row[f_list.index("Evidence")]
        pmid = "29351928"
        if ac in ac2canon:
            canon = ac2canon[ac]
            for i in xrange(0, len(s_list)):
                s, compo = s_list[i], g_list[i]
                if s == "":
                    continue
                aa_one, aa_pos = s[0], int(s[1:])
                if aa_one not in aa_one2three:
                    continue
                aa_three = aa_one2three[aa_one]
                if aa_pos == 0 or aa_pos >= len(seq_hash[canon]):
                    continue
                if aa_one != seq_hash[canon][aa_pos-1]:
                    continue
                newrow_list = [
                    [ac2canon[ac],"o-linked","G70994MS", aa_three, str(aa_pos), evidence, compo],
                    [ac2canon[ac],"o-linked","G70994MS", aa_three, str(aa_pos), pmid, compo]
                ]
                for newrow in newrow_list:
                    newrowstr = json.dumps(newrow)
                    if newrowstr not in seen:
                        print "\"%s\"" % ("\",\"".join(newrow))
                        seen[newrowstr] = True
    return

def load_glycosylation_type_one(species):


    data_frame = {}
    in_file = path_obj["downloads"] + "unicarbkb/%s_motif_current.txt" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    glytoucanac2glycosylationtype = {}
    for row in data_frame["data"]:
        uckb_id = row[f_list.index("uckb_id")].strip()
        glytoucan_ac = row[f_list.index("motif_ac")].strip()
        motif_label = row[f_list.index("motif_label")].strip()
        if glytoucan_ac not in glytoucanac2glycosylationtype:
            glytoucanac2glycosylationtype[glytoucan_ac] = []
        if motif_label.lower().find("n-glycan") != -1:
            if "n-linked" not in glytoucanac2glycosylationtype[glytoucan_ac]:
                glytoucanac2glycosylationtype[glytoucan_ac].append("n-linked")
        if motif_label.lower().find("o-glycan") != -1:
            if "o-linked" not in glytoucanac2glycosylationtype[glytoucan_ac]:
                glytoucanac2glycosylationtype[glytoucan_ac].append("o-linked")
    
    return glytoucanac2glycosylationtype




def load_glycosylation_type_two():

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/classification.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    glytoucanac2glycosylationtype = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")].strip()
        gly_type = row[f_list.index("Type")].strip().lower()
        gly_subtype = row[f_list.index("Subtype")].strip()
        if glytoucan_ac not in glytoucanac2glycosylationtype:
            glytoucanac2glycosylationtype[glytoucan_ac] = []
        if gly_type not in glytoucanac2glycosylationtype[glytoucan_ac]:
            if gly_type == "n-linked":
                glytoucanac2glycosylationtype[glytoucan_ac].append(gly_type)
            if gly_type == "o-linked":
                glytoucanac2glycosylationtype[glytoucan_ac].append(gly_type)

    return glytoucanac2glycosylationtype





def extract_glycosylation_sites_unicarbkb_ds(species):


    #glytoucanac2glycosylationtype = load_glycosylation_type_one(species)
    glytoucanac2glycosylationtype = load_glycosylation_type_two()


    #We do stringent ac2canon and insist that isoform_ac is the same as canon_ac
    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]] 
        for isoform in isoform_list:
            isoform_ac = isoform.split("-")[0].strip()
            canon_ac = canon.split("-")[0].strip()
            if canon_ac == isoform_ac:
                ac2canon[isoform_ac] = canon


    glycan_list = []
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/glycan_properties.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_acc")]
        if ac not in glycan_list:
            glycan_list.append(ac)

    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_allsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())


    n_linked_aalist = []
    o_linked_aalist = []
    
    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_three2one = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        gly_type = row[f_list.index("gly_type")] 
        if gly_type != "":
            aa_three2one[three] = one
            if gly_type == "n-linked":
                n_linked_aalist.append(three)
            elif gly_type == "o-linked":
                o_linked_aalist.append(three)


    log_file_one = path_obj["logs"] + "%s_proteoform_glycosylation_sites_unicarbkb.1.log" % (species)
    FL1 = open(log_file_one, "w")
    log_file_two = path_obj["logs"] + "%s_proteoform_glycosylation_sites_unicarbkb.2.log" % (species)
    FL2 = open(log_file_two, "w")

    log_file_three = path_obj["intermediate"] + "%s_proteoform_glycosylation_sites_unicarbkb.csv" % (species)
    FL3 = open(log_file_three, "w")


    newrow = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","evidence",
            "unicarbkb_id","saccharide","amino_acid","glycosylation_type", "curation_notes", 
            "additonal_notes","composition"]
    print "\"%s\"" % ("\",\"".join(newrow))
    FL3.write("\"%s\"\n" % ("\",\"".join(newrow + ["qc_tag"])))


    data_frame = {}
    #in_file = path_obj["downloads"] + "unicarbkb/%s_glycosylation_current.csv" % (species)
    #in_file = "/data/projects/glygen/downloads/unicarbkb/human29112019.csv"
    in_file = "/data/projects/glygen/downloads/unicarbkb/human17122019.csv"



    
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    n, n1, n2 = 0, 0, 0

    for row in data_frame["data"]:
        uckb_id = row[f_list.index("Id")]
        uniprotkb_ac = row[f_list.index("Protein")].split("/")[-1]
        glytoucan_ac = row[f_list.index("Toucan")] if "Toucan" in f_list else ""
        pmid = row[f_list.index("Pmid")]
        pos = row[f_list.index("Position")].split("^^")[0]
        pos = int(pos) if pos != "" else 0
        amino_acid = row[f_list.index("TypeAminoAcid")]
        amino_acid = amino_acid[0].upper() +  amino_acid[1:].lower()
        notes = row[f_list.index("Notes")]
        additional_notes = row[f_list.index("AdditionalNotes")] if "AdditionalNotes" in f_list else ""
        
        composition = ""
        flag_dict = {
                "no_canon":False, "bad_pos":False, "aa_mismatch":False, "invalid_aa":False,
                "glycan_without_glytype":False,"bad_glytype":False,"no_pmid":False, 
                "no_glycan_invalid_aa":False, "glycan_without_n_or_o_glytype": False
        }
        if uniprotkb_ac not in ac2canon:
            flag_dict["no_canon"] = True
            canon = uniprotkb_ac
        else:
            canon = ac2canon[uniprotkb_ac]
            if amino_acid in aa_three2one:
                if pos == 0 or pos >= len(seq_hash[canon]):
                    flag_dict["bad_pos"] = True
                elif aa_three2one[amino_acid] != seq_hash[canon][pos-1]:
                    flag_dict["aa_mismatch"] = True
            else: #### amino acid is not n or o linked
                flag_dict["invalid_aa"] = True
        
        glycosylation_type = "xxx: should change by under one of the conditions"
        if glytoucan_ac == "":
            if amino_acid in n_linked_aalist:
                glycosylation_type = "n-linked"
            elif amino_acid in o_linked_aalist:
                glycosylation_type = "o-linked"
            else:
                glycosylation_type = ""
                flag_dict["no_glycan_invalid_aa"] = True
        elif glytoucan_ac not in glytoucanac2glycosylationtype:
            glycosylation_type = ""
            flag_dict["glycan_without_glytype"] = True
        elif glytoucanac2glycosylationtype[glytoucan_ac] == []:
            glycosylation_type = ""
            #flag_dict["glycan_without_n_or_o_glytype"] = True
        else: #resume normal function if we have reported glycosylation information 
            tmp_list = []
            if glytoucan_ac in glytoucanac2glycosylationtype:
                if "n-linked" in glytoucanac2glycosylationtype[glytoucan_ac]:
                    tmp_list.append("n-linked")
                if "o-linked" in glytoucanac2glycosylationtype[glytoucan_ac]:
                    tmp_list.append("o-linked")
                if tmp_list == ["n-linked"] and amino_acid not in n_linked_aalist:
                    flag_dict["n_glycan_aa_mismatch"] = True
                if tmp_list == ["o-linked"] and amino_acid not in o_linked_aalist:
                    flag_dict["o_glycan_aa_mismatch"] = True
            glycosylation_type = ";".join(tmp_list)
            if glycosylation_type == "n-linked;o-linked":
                flag_dict["bad_glytype"] = True
        
        if pmid in ["", "0"]:
            flag_dict["no_pmid"] = True
        
        if glytoucan_ac == "" and uckb_id[0:5] == "comp_":
            composition = uckb_id
        if uckb_id[0:5] == "comp_":
            uckb_id = ""
        if glytoucan_ac != "" and glytoucan_ac not in glycan_list:
            FL2.write("\"%s\"\n" % (glytoucan_ac))
    
        newrow = [canon,str(pos),pmid, uckb_id, glytoucan_ac,amino_acid,glycosylation_type, notes,
            additional_notes, composition]

        n += 1
        flag_list = []
        for flag in flag_dict:
            if flag_dict[flag] == True:
                flag_list.append(flag)
        if flag_list != []:
            newrow.append(";".join(flag_list))
            FL1.write("\"%s\"\n" % ("\",\"".join(newrow)))
            FL3.write("\"%s\"\n" % ("\",\"".join(newrow)))
            n2 += 1
        else:
            n1 += 1
            FL3.write("\"%s\"\n" % ("\",\"".join(newrow + ["validation_passed"])))
            print "\"%s\"" % ("\",\"".join(newrow))
    FL1.close()
    FL2.close()
    FL3.close()

    return



def extract_phosphorylation_sites_uniprotkb_ds(species):

    data_grid = { "ac2canon":{}, "isoformlist":{},  "phosphorylation":[], "genename":{}}
    sparqlutil.load_phosphorylation_sites(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_genename(data_grid, species)
    
    row = ["uniprotkb_canonical_ac", "gene_symbol","phosphorylated_residue","amino_acid",
            "phosphorylation_site_uniprotkb", "protein_kinase","data_source", "evidence","eco_id"]
    print "\"%s\""  % ("\",\"".join(row))
    for row in data_grid["phosphorylation"]:
        if row[0] in data_grid["ac2canon"]:
            ac = row[0]
            gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
            newrow = [data_grid["ac2canon"][ac], gene_name] + row[1:]
            print "\"%s\""  % ("\",\"".join(newrow))
    return



def extract_glycosylation_sites_uniprotkb_ds(species):

    data_grid = { "ac2canon":{}, "isoformlist":{},  "glycosylation":[]}
    sparqlutil.load_glycosylation_sites(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)

    row = ["uniprotkb_canonical_ac","glycosylation_type","saccharide","amino_acid",
            "glycosylation_site_uniprotkb","uniprotkb_glycosylation_annotation_comment", 
            "data_source","evidence","eco_id","uniprotkb_ftid"]
    print "\"%s\""  % ("\",\"".join(row))
    for row in data_grid["glycosylation"]:
        if row[0] in data_grid["ac2canon"]:
            row[0] = data_grid["ac2canon"][row[0]]
            
            #Ignore glycation cases
            if row[1].lower() == "n-linked" and row[3] != "Asn":
                continue
            
            print "\"%s\""  % ("\",\"".join(row))

    return


def extract_glycosylation_sites_pdb_ds(species):

    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_one2three = {}
    aa_three2one = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        aa_one2three[one] = three
        aa_three2one[three] = one

    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon

    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_canonicalsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())

    data_frame = {}
    in_file = "downloads/pdb/current/pdb_chain_uniprot.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    chain_dict = {}
    for row in data_frame["data"]:
        pdb_id = row[f_list.index("PDB")]
        pdb_chain = row[f_list.index("CHAIN")]
        ac = row[f_list.index("SP_PRIMARY")]
        res_beg = row[f_list.index("RES_BEG")]
        sp_beg =row[f_list.index("SP_BEG")]
        combo = "%s:%s" % (pdb_id.lower(), pdb_chain.lower())
        o = {"ac":ac, "resbeg":res_beg, "spbeg":sp_beg}
        if combo not in chain_dict:
            chain_dict[combo] = []
        chain_dict[combo].append(o)

    newrow = ["uniprotkb_canonical_ac","glycosylation_type","saccharide","amino_acid",
            "glycosylation_site_uniprotkb","evidence",
            "pdb_id","amino_acid_chain","glycosylation_site_pdb"
    ]
    print "\"%s\""  % ("\",\"".join(newrow))

    for in_file in glob.glob(path_obj["downloads"] + "pdb/current/*-linked-*-details"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, "\t")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            aathree_in_pdb = in_file.split("/")[-1].split("-")[-2]
            glycosylation_type = in_file.split("/")[-1].split("-")[0] + "-linked"
            pdb_id = row[f_list.index("PDB-ID")]
            pdb_chain = row[f_list.index("PDB-Chain")]
            amino_acid = row[f_list.index("Amino-Acid")]
            pdb_pos = int(row[f_list.index("PDB-Position")])
            saccharide_name = row[f_list.index("Carb")]
            combo = "%s:%s" % (pdb_id.lower(), pdb_chain.lower())
            if combo not in chain_dict:
                continue
            if len(chain_dict[combo]) > 1:
                continue
            o = chain_dict[combo][0]
            ac, sp_beg, res_beg = o["ac"], int(o["spbeg"]), int(o["resbeg"])
            canon_pos = pdb_pos + (sp_beg - res_beg)
            if ac not in ac2canon:
                continue
            canon = ac2canon[ac]
            if canon_pos < 0 or canon_pos > len(seq_hash[canon]):
                continue
            aaone_in_canon = seq_hash[canon][canon_pos-1]
            aathree_in_pdb = aathree_in_pdb[0].upper() + aathree_in_pdb[1:].lower()
            aaone_in_pdb = aa_three2one[aathree_in_pdb]
            if aaone_in_canon != aaone_in_pdb:
                continue
            newrow = [canon,glycosylation_type,saccharide_name,aathree_in_pdb,str(canon_pos),"",
                    pdb_id,pdb_chain,str(pdb_pos)]
            print "\"%s\""  % ("\",\"".join(newrow))


    return


def extract_glycosylation_sites_pdb_ds_old(species):


    pdbid_list = []
    for f in glob.glob(path_obj["downloads"] + "pdb/current/pdb_format/*.pdb"):
        pdb_id = f.split("/")[-1].replace(".pdb", "")
        if pdb_id not in pdbid_list:
            pdbid_list.append(pdb_id)


    log_file = path_obj["logs"] + "proteoform_glycosylation_sites_pdb.log"
    out_file = path_obj["unreviewed"] + "proteoform_glycosylation_sites_pdb.csv"
    
    distmin = 1.35
    distmax = 1.50 
    header_list = ["uniprotkb_acc","pdbid_list","amino_acid","glycosylation_type","amino_acid_chain","glycosylation_site","glycosylation_site_uniprotkb"]
    row_list = []
    row_list.append(header_list)

    FL = open(log_file, "w")
    for id in pdbid_list:
        flag_dbref = 0
        flag_dbref1 = 0
        flag_modres = 0
        flag_link = 0
        rowcount=0
        acc_glycosylation = []
        dbref = {}
        in_file = path_obj["downloads"] + "pdb/current/pdb_format/%s.pdb" % (id)
        FR = open(in_file, "r") 
        for row in FR:
            rowcount+=1
            line = " ".join(row.split()).split(" ")
            if "DBREF" in line[0]:
                cnd_list = ["UNP" in line[5], len(line)==10, line[-2].replace("-","").isdigit()]
                cnd_list.append(line[3].replace("-","").isdigit())
                if False not in cnd_list:
                    flag_dbref = 1
                    if line[2] not in dbref:
                        dbref[line[2]]=[]
                    dbref[line[2]]=[line[6],int(line[-2])-int(line[3])]
                elif "UNP" in line[5] and len(line)==7 and line[3].replace("-","").isdigit():
                    if line[2] not in dbref:
                        flag_dbref1 = 1
                        temp_chain = line[2]
                        glycosylation_site_pdb = line[3]
                elif len(line)==6 and line[2]==temp_chain:
                    #print line
                    dbref[line[2]]=[line[3],int(line[-2])-int(glycosylation_site_pdb)]
                    flag_dbref1 = 0
                else:
                    continue
            if "MODRES" in line[0] and "GLYCOSYLATION" in line[-2] and "SITE" in line[-1]:
                flag_modres = 1
                if "ASN" in line and line[line.index("ASN")+2].replace("-","").isdigit() and line[line.index("ASN")+1] in dbref:
                    chain = line[line.index("ASN")+1]
                    glycosylation_site = int(line[line.index("ASN")+2])
                    if [dbref[chain][0],glycosylation_site] not in acc_glycosylation:
                        acc_glycosylation.append([dbref[chain][0],glycosylation_site])
                        tmprow = [dbref[chain][0],id,"ASN","N-Linked",chain,str(glycosylation_site),
                                str(glycosylation_site+dbref[chain][1])]
                        row_list.append(tmprow)
                elif "THR" in line and line[line.index("THR")+2].replace("-","").isdigit() and line[line.index("THR")+1] in dbref:
                    chain = line[line.index("THR")+1]
                    glycosylation_site = int(line[line.index("THR")+2])
                    if [dbref[chain][0],glycosylation_site] not in acc_glycosylation:
                        acc_glycosylation.append([dbref[chain][0],glycosylation_site])
                        tmprow = [dbref[chain][0],id,"THR","O-Linked",chain,str(glycosylation_site),
                                str(glycosylation_site+dbref[chain][1])]
                        row_list.append(tmprow)
                elif "SER" in line and line[line.index("SER")+2].replace("-","").isdigit() and line[line.index("SER")+1] in dbref:
                    chain = line[line.index("SER")+1]
                    glycosylation_site = int(line[line.index("SER")+2])
                    if [dbref[chain][0],glycosylation_site] not in acc_glycosylation:
                        acc_glycosylation.append([dbref[chain][0],glycosylation_site])
                        tmprow = [dbref[chain][0],id,"SER","O-Linked",chain,str(glycosylation_site),
                                str(glycosylation_site+dbref[chain][1])]
                        row_list.append(tmprow)
                else:
                    flag_modres = 0
                    continue
            if "LINK" in line[0] and flag_modres ==0: 
                if row[16] != " ":
                    row = row[:16] + ' ' + row[17:]
                if row[22].replace("-","").isdigit():
                    row = row[:22] + ' ' + row[22:]
                if row[46] != " ":
                    row = row[:46] + ' ' + row[47:]
                if row[52].replace("-","").isdigit():
                    row = row[:52] + ' ' + row[52:]
                line = " ".join(row.split()).split(" ")
                if "ASN" in line and line[line.index("ASN")+2].replace("-","").isdigit() and line[line.index("ASN")+1] in dbref:
                    dist = float(line[-1])
                    if distmin<=dist and dist<=distmax:
                        flag_link = 1
                        chain = line[line.index("ASN")+1]
                        glycosylation_site = int(line[line.index("ASN")+2])
                        if [dbref[chain][0],glycosylation_site] not in acc_glycosylation:
                            acc_glycosylation.append([dbref[chain][0],glycosylation_site])
                            tmprow = [dbref[chain][0],id,"ASN","N-Linked",
                                    chain,str(glycosylation_site),str(glycosylation_site+dbref[chain][1])]
                            row_list.append(tmprow)
                    else: 
                        FL.write(id+"\t"+ "Distance is not located in the range between 1.35 and 1.5."+"\t"+
                                " ".join(row.split())+"\n")
                elif "THR" in line and line[line.index("THR")+2].replace("-","").isdigit() and line[line.index("THR")+1] in dbref:
                    dist = float(line[-1])
                    if distmin<=dist and dist<=distmax:
                        flag_link = 1
                        chain = line[line.index("THR")+1]
                        glycosylation_site = int(line[line.index("THR")+2])
                        if [dbref[chain][0],glycosylation_site] not in acc_glycosylation:
                            acc_glycosylation.append([dbref[chain][0],glycosylation_site])
                            tmprow = [dbref[chain][0],id,"THR","O-Linked",
                                    chain,str(glycosylation_site),str(glycosylation_site+dbref[chain][1])]
                            row_list.append(tmprow)
                    else:
                        FL.write(id+"\t"+ "Distance is not located in the range between 1.35 and 1.5."+"\t" + 
                                " ".join(row.split())+"\n")
                elif "SER" in line and line[line.index("SER")+2].replace("-","").isdigit() and line[line.index("SER")+1] in dbref:
                    dist = float(line[-1])
                    if distmin<=dist and dist<=distmax:
                        chain = line[line.index("SER")+1]
                        glycosylation_site = int(line[line.index("SER")+2])
                        if [dbref[chain][0],glycosylation_site] not in acc_glycosylation:
                            acc_glycosylation.append([dbref[chain][0],glycosylation_site])
                            tmprow = [dbref[chain][0],id,"SER","O-Linked",
                                    chain,str(glycosylation_site),str(glycosylation_site+dbref[chain][1])]
                            row_list.append(tmprow)
                    else:
                        FL.write(id+"\t"+ "Distance is not located in the range between 1.35 and 1.5."+"\t"+ 
                                " ".join(row.split())+"\n")
                else:
                    flag_link = 0
                    continue
        if flag_dbref == 0:
            FL.write(id+"\t"+ "There is not uniprot id in DBREF record."+"\n")
            continue
        if flag_modres == 0 and flag_link == 0:
            FL.write(id+"\t"+"There are not any records about modres and link."+"\n")
            continue
        FR.close()
    FL.close()


    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_canonicalsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())
        
    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon

    newrow = ["uniprotkb_canonical_ac","pdb_id","amino_acid","glycosylation_type",
        "amino_acid_chain","glycosylation_site_pdb","glycosylation_site_uniprotkb","amino_acid_uniprotkb"
    ]
    print "\"%s\"" % ("\",\"".join(newrow))

    for row in row_list:
        if row[0] not in ac2canon:
            continue
        canon = ac2canon[row[0]]
        pos = int(row[6])
        if pos < 0 or pos > len(seq_hash[canon]):
            continue
        aa_in_canon = seq_hash[canon][pos-1]
        aa_in_pdb = row[2]
        newrow = row
        newrow[0] = canon
        if aa_in_pdb == "ASN" and aa_in_canon == "N":
            newrow += ["Asn"]
        elif aa_in_pdb == "SER" and aa_in_canon == "S":
            newrow += ["Ser"]
        elif aa_in_pdb == "Thr" and aa_in_canon == "T":
            newrow += ["Thr"]
        else:
            continue
        print "\"%s\"" % ("\",\"".join(newrow))

    return


###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")
    parser.add_option("-d","--dataset",action="store",dest="dataset",help="[glycosylation_sites_unicarbkb, glycosylation_sites_pdb]")
    (options,args) = parser.parse_args()
    for file in ([options.species, options.dataset]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    global config_obj
    global sparql
    global graph_uri
    global prefixes
    global data_grid
    global species_obj


    species = options.species
    dataset = options.dataset

    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]
    
    global path_obj
    path_obj = config_obj["pathinfo"]


    data_grid = {}
    if dataset == "glycosylation_sites_unicarbkb":
        extract_glycosylation_sites_unicarbkb_ds(species)
    elif dataset == "glycosylation_sites_uniprotkb":
        extract_glycosylation_sites_uniprotkb_ds(species)
    elif dataset == "glycosylation_sites_pdb":
        extract_glycosylation_sites_pdb_ds(species)
    elif dataset == "glycosylation_sites_harvard":
        extract_glycosylation_sites_harvard_ds(species)
    elif dataset == "glycosylation_sites_literature":
        extract_glycosylation_sites_literature_ds(species)
    elif dataset == "phosphorylation_sites_uniprotkb":
        extract_phosphorylation_sites_uniprotkb_ds(species)
    elif dataset == "citations_unicarbkb":
        extract_citations_unicarbkb_ds(species)
    elif dataset == "glycosylation_sites_tyr_o_linked":
        extract_glycosylation_sites_tyr_o_linked_ds(species)
    elif dataset == "glycosylation_sites_glyconnect":
        extract_glycosylation_sites_glyconnect_ds(species)







if __name__ == '__main__':
        main()

