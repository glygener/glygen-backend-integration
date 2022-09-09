import os,sys
import json
import csv
import time

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
#import requests
import gzip


import libgly


__version__="1.0"
__status__ = "Dev"


def extract_top_authors_ds():

    file_list = glob.glob("unreviewed/glycan_citations_*.csv")
    seen = {}
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            glytoucan_ac = row[f_list.index("glytoucan_ac")]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            author_list = row[f_list.index("authors")].split(",")
            for author in author_list:
                author = author.strip()
                if author not in seen:
                    seen[author] = {}
                if xref_id not in seen[author]:
                    seen[author][xref_id] = {}
                seen[author][xref_id][glytoucan_ac] = True

    row_dict = {}
    for author in seen:
        pmid_list = []
        glycan_list = []
        for xref_id in seen[author]:
            pmid_list.append(xref_id)
            for glytoucan_ac in seen[author][xref_id]:
                glycan_list.append(glytoucan_ac)
        newrow = [author, str(len(glycan_list)), "|".join(glycan_list), "|".join(pmid_list)] 
        n = len(glycan_list)
        if n not in row_dict:
            row_dict[n] = []
        row_dict[n].append(newrow)

    newrow = ["author_name", "glycan_count", "glytoucan_ac_list"    ,"pmid_list"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for n in sorted(row_dict, reverse=True):
        for newrow in row_dict[n]:
            print "\"%s\"" % ("\",\"".join(newrow))


    return



def extract_glygen_uniprotkb_protvista_mapping_ds():
    
    glycan_list = load_glycan_masterlist()

    src_list_one = ["unicarbkb", "glyconnect", "gptwiki", "harvard", "o_glcnac_mcw"]
    src_list_two = ["uniprotkb"]
    f_list_one = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","amino_acid","saccharide","glycosylation_type","xref_key","xref_id","src_xref_key","src_xref_id"]
    f_list_two = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","xref_key","xref_id","src_xref_key","src_xref_id"]

    obj_dict_one = {}
    for src in src_list_one:
        file_list = glob.glob("unreviewed/*_proteoform_glycosylation_sites_%s.csv" % (src))
        for in_file in file_list:
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            f_list = data_frame["fields"]
            species = in_file.split("/")[-1].split("_")[0]
            tax_id = species_obj[species]["tax_id"]
            for row in data_frame["data"]:
                ac = row[f_list.index("uniprotkb_canonical_ac")].split("-")[0]
                obj = {"tax_id":tax_id, "uniprotkb_ac":ac}
                for f in f_list_one:
                    obj[f]= row[f_list.index(f)]
                if obj["saccharide"] not in glycan_list:
                    continue
                canon,pos = obj["uniprotkb_canonical_ac"],obj["glycosylation_site_uniprotkb"]
                site_id = "%s %s" % (canon,pos)
                if site_id not in obj_dict_one:
                    obj_dict_one[site_id] = []
                obj_dict_one[site_id].append(obj)
    
    obj_dict_two = {}
    for src in src_list_two:
        file_list = glob.glob("unreviewed/*_proteoform_glycosylation_sites_%s.csv" % (src))
        for in_file in file_list:
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            f_list = data_frame["fields"]
            species = in_file.split("/")[-1].split("_")[0]
            tax_id = species_obj[species]["tax_id"]
            for row in data_frame["data"]:
                ac = row[f_list.index("uniprotkb_canonical_ac")].split("-")[0]
                obj = {"tax_id":tax_id, "uniprotkb_ac":ac}
                for f in f_list_two:
                    obj[f]= row[f_list.index(f)]
                canon,pos = obj["uniprotkb_canonical_ac"],obj["glycosylation_site_uniprotkb"]
                site_id = "%s %s" % (canon,pos)
                if site_id not in obj_dict_two:
                    obj_dict_two[site_id] = []
                obj_dict_two[site_id].append(obj)


    gtc2chebi = {}
    in_file = "unreviewed/glycan_xref_chebi.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")]
        gtc2chebi[glytoucan_ac] = row[f_list.index("xref_id")]
    



    newrow =  ["uniprotkb_ac"] + f_list_one[:-4] + ["chebi_id", "tax_id","xref_key","xref_id",
            "src_xref_key","src_xref_id"]
    row_str = "\t".join(newrow)
    print row_str

    seen = {}
    for site_id in obj_dict_two:
        if site_id not in obj_dict_one:
            continue
        for obj in obj_dict_one[site_id]:
            if obj["saccharide"] not in gtc2chebi:
                continue
            obj["chebi_id"] = gtc2chebi[obj["saccharide"]]
            pair_list = []
            combo_id = "%s %s %s %s" % (obj["xref_key"],obj["xref_id"],obj["src_xref_key"],obj["src_xref_id"])
            if combo_id not in pair_list:
                pair_list.append(combo_id)
            for o in obj_dict_two[site_id]:
                combo_id = "%s %s %s %s" % (o["xref_key"],o["xref_id"],o["src_xref_key"],o["src_xref_id"])
                if combo_id not in pair_list:
                    pair_list.append(combo_id)
            row = []
            for f in ["uniprotkb_ac"] + f_list_one[:-4] + ["chebi_id", "tax_id"]:
                row.append(str(obj[f]))
            for combo_id in pair_list:
                xref_key, xref_id, src_xref_key, src_xref_id = combo_id.split(" ")
                if xref_key in ["protein_xref_uniprotkb_gly"]:
                    continue
                xref_id = xref_id.replace("GLYDS", "GLY_")
                src_xref_id = src_xref_id.replace("GLYDS", "GLY_")
                newrow = row + [xref_key, xref_id, src_xref_key, src_xref_id]
                row_str = "\t".join(newrow)
                if row_str not in seen:
                    print row_str
                seen[row_str] = True


    #Map glytoucan_ac to chebi_id using file 8.
    #Keep rows that have chebi_id
    #*Format of the output file:* (should be tsv- as requested by UniProtKB)
    #"uniprotkb_ac","uniprotkb_canonical_ac","glycosylation_site_uniprotkb","amino_acid","saccharide","glycosylation_type","xref_key","xref_id","tax_id","chebi_id"

    return



def load_motif_masterlist():

    motif_list = []
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        is_motif = row[f_list.index("is_motif")]
        if is_motif == "True":
            motif_list.append(ac)

    return motif_list



def load_glycan_masterlist():

    glycan_list = []
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        if ac not in glycan_list:
            glycan_list.append(ac)

    return glycan_list




def get_blacklisted_pmids():

    black_list = []
    in_file = "compiled/%s_protein_blacklisted_pmids.csv" 
    if os.path.isfile(in_file) == True:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        for row in data_frame["data"]:
            black_list.append(row[0])
        black_list = sorted(set(black_list))

    return black_list




def extract_sequences_ds(seq_format):


    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/monocomp.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    ac2comp = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("accession")]
        comp = ""
        for f in f_list[1:-1]:
            n = row[f_list.index(f)]
            if n != "0":
                comp += "%s(%s)" % (f, n)
        ac2comp[glytoucan_ac] = comp


    newrow = ["glytoucan_ac","sequence_%s" % (seq_format)]
    print "\"%s\"" % ("\",\"".join(newrow))

    seen_list = []
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        mass = row[f_list.index("glycan_mass")]
        if ac not in seen_list:
            seen_list.append(ac)
            if seq_format == "byonic":
                seq = "%s %s %s" % (ac2comp[ac], "%", mass) if ac in ac2comp else ""
                if seq.find("Xxx") == -1:
                    newrow = [ac, seq]
                    print "\"%s\"" % ("\",\"".join(newrow))
            else:
                in_file = path_obj["downloads"] + "glytoucan/current/export/%s/%s.txt" % (seq_format, ac)
                if os.path.isfile(in_file) == True:
                    with open(in_file, "r") as FR:
                        seq = ""
                        for line in FR:
                            seq += " " + line.strip()
                        newrow = [ac, seq.strip()]
                        print "\"%s\"" % ("\",\"".join(newrow))

    return



    
def extract_monosaccharide_composition_ds():

    glycan_list = load_glycan_masterlist()

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/monocomp.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    newrow = f_list
    newrow[0] = "glytoucan_ac"
    print "\"%s\"" % ("\",\"".join(newrow))

    seen = {}
    for row in data_frame["data"]:
        newrow = row
        if newrow[0] in glycan_list:
            print "\"%s\"" % ("\",\"".join(newrow))

    return


def extract_monosaccharide_composition_advanced_ds():

    glycan_list = load_glycan_masterlist()

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/monocounts.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    count_dict = {}
    seen_mono = {}
    for row in data_frame["data"]:
        newrow = row
        ac = row[f_list.index("accession")]
        monosaccharide = row[f_list.index("monosaccharide")]
        #monosaccharide = "Xxx" if monosaccharide == "Pent" else monosaccharide
        count = row[f_list.index("count")].replace("+", "")
        if ac in glycan_list:
            if ac not in count_dict:
                count_dict[ac] = {}
            count_dict[ac][monosaccharide] = int(count)
            if monosaccharide != "*":
                seen_mono[monosaccharide] = True


    mono_list = sorted(seen_mono.keys()) + ["Count"]
    newrow = ["glytoucan_ac"] + mono_list
    print "\"%s\"" % ("\",\"".join(newrow))
    for ac in count_dict:
        newrow = [ac]
        for monosaccharide in mono_list:
            n = count_dict[ac][monosaccharide] if monosaccharide in count_dict[ac] else 0
            newrow.append(str(n))
        print "\"%s\"" % ("\",\"".join(newrow))

    return



def extract_glytoucan_accession_history_ds():
    
    glycan_list = load_glycan_masterlist()
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/glygen_retired_accessions.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    newrow = ["glytoucan_ac_old","glytoucan_ac_current"]
    print "\"%s\"" % ("\",\"".join(newrow))

    for row in data_frame["data"]:
        glytoucan_ac_old = row[f_list.index("accession")]
        glytoucan_ac_current = row[f_list.index("replacewith")]
        glytoucan_ac_current = glytoucan_ac_current if glytoucan_ac_current != "-" else "" 
        newrow = [glytoucan_ac_old, glytoucan_ac_current]
        print "\"%s\"" % ("\",\"".join(newrow))

    return 


def extract_motif_ds():


    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/allmotifs.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    
    motif_dict = {}
    for row in data_frame["data"]:
        motif_ac = row[f_list.index("MotifAccession")]
        prop = row[f_list.index("Property")]
        value = row[f_list.index("Value")]
        if motif_ac not in motif_dict:
            motif_dict[motif_ac] = {}
        if prop not in motif_dict[motif_ac]:
            motif_dict[motif_ac][prop] = value
        else:
            motif_dict[motif_ac][prop] += "; " + value


    glycan_list = load_glycan_masterlist()
    newrow = ["glytoucan_ac","motif_ac","motif_name","alignment","collection_name",
        "collection_accession", "alternative_name", "aglycon", "motif_ac_xref",
        "reducing_end", "pmid", "keyword"
    ]
    print "\"%s\"" % ("\",\"".join(newrow))
    prop_list = [
        "PreferredName", "Alignment", "Collection", "Accession","AlternativeName",
        "Aglycon", "MotifGlyTouCan","ReducingEnd", "PMID", "Keyword"
    ]

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/allmotifaligns.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    seen = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        motif_ac = row[f_list.index("MotifAccession")]
        motif_label = row[f_list.index("Label")]
        #is_reducing_end = row[f_list.index("IsReducingEnd")]
        alignment = row[f_list.index("Alignment")]

        combo_id = "%s %s" % (ac, motif_ac)
    
        p = "MotifGlyTouCan"
        motif_ac_xref = motif_dict[motif_ac][p] if p in motif_dict[motif_ac] else ""
        if ac in glycan_list and motif_ac_xref in glycan_list:
            if combo_id not in seen:
                seen[combo_id] = True
                newrow = [ac, motif_ac]
                for p in prop_list:
                    val = motif_dict[motif_ac][p] if p in motif_dict[motif_ac] else ""
                    newrow.append(val)
                print "\"%s\"" % ("\",\"".join(newrow))




    return



def extract_classification_ds():

    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","glycan_type","glycan_type_source", "glycan_type_source_id", 
            "glycan_subtype", "glycan_subtype_source", "glycan_subtype_source_id"]
    print "\"%s\"" % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/classification.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    type_dict = {}
    subtype_dict = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        glycan_class = row[f_list.index("Classification")]
        level = row[f_list.index("Level")]
        if glytoucan_ac not in type_dict:
            type_dict[glytoucan_ac] = {}
        if level == "GlycanType":
            glycan_type = glycan_class
            source = row[f_list.index("Source")]
            source_id = row[f_list.index("SourceID")]
            o = {"name":glycan_type, "source":source, "source_id":source_id}
            if glycan_type not in type_dict[glytoucan_ac]:
                type_dict[glytoucan_ac][glycan_type] = {}
            type_dict[glytoucan_ac][glycan_type] = o
        elif level == "GlycanSubtype":
            source = row[f_list.index("Source")]
            source_id = row[f_list.index("SourceID")]
            glycan_subtype = glycan_class
            glycan_type = glycan_class.split(" ")[0]
            if glytoucan_ac not in subtype_dict:
                subtype_dict[glytoucan_ac] = {}
            o = {"name":glycan_subtype, "source":source, "source_id":source_id}
            if glycan_type not in subtype_dict[glytoucan_ac]:
                subtype_dict[glytoucan_ac][glycan_type] = {}
            subtype_dict[glytoucan_ac][glycan_type][glycan_subtype] = o
    

    for glytoucan_ac in type_dict:
        if glytoucan_ac in glycan_list:
            for glycan_type in type_dict[glytoucan_ac]:
                o = type_dict[glytoucan_ac][glycan_type]
                t_name, t_src, t_src_id = o["name"], o["source"], o["source_id"]
                if glytoucan_ac in subtype_dict:
                    if glycan_type in subtype_dict[glytoucan_ac]:
                        for glycan_subtype in subtype_dict[glytoucan_ac][glycan_type]:
                            oo = subtype_dict[glytoucan_ac][glycan_type][glycan_subtype]
                            s_name, s_src, s_src_id = oo["name"], oo["source"], oo["source_id"]
                            s_name =  " ".join(s_name.split(" ")[1:]).strip()
                            s_name = s_name[0].upper() + s_name[1:]
                            newrow = [glytoucan_ac,t_name,t_src,t_src_id,s_name,s_src,s_src_id]
                            print "\"%s\"" % ("\",\"".join(newrow))
                else:
                    s_name, s_src, s_src_id = "", "", ""
                    newrow = [glytoucan_ac,t_name,t_src,t_src_id,s_name,s_src,s_src_id]
                    print "\"%s\"" % ("\",\"".join(newrow))
                
        

    return

def extract_citations_glycomotif_ds():
    black_list = get_blacklisted_pmids()
    glycan_list = load_glycan_masterlist()

    in_file = path_obj["unreviewed"] + "/glycan_motif.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    
    FL = open(path_obj["logs"] + "/glycan_citations_glycomotif.log", "w")
    newrow = ["glytoucan_ac","title","journal_name","publication_date", "authors",
            "xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("motif_ac_xref")]
        motif_ac = row[f_list.index("motif_ac")]
        for xref_id in row[f_list.index("pmid")].strip().split(" "):
            xref_id = xref_id.replace(";", "")
            if xref_id == "":
                continue
            xref_key = "glycan_xref_pubmed"
            src_xref_id = motif_ac
            src_xref_key = "glycan_xref_motif"
            combo_id = "%s %s" % (glytoucan_ac, xref_id)
            if glytoucan_ac not in glycan_list:
                FL.write("%s\n" %(glytoucan_ac))
            citerow = libgly.get_citation(xref_id, path_obj["downloads"] + "ncbi/medline/")

            if citerow != [] and combo_id not in seen:
                newrow = [glytoucan_ac] + citerow + [xref_key,xref_id, src_xref_key,src_xref_id]
                print "\"%s\"" % ("\",\"".join(newrow))
                seen[combo_id] = True

    FL.close()




def extract_citations_ncfg_ds():
    black_list = get_blacklisted_pmids()
    glycan_list = load_glycan_masterlist()

    in_file = path_obj["downloads"] + "ncfg/current/glycan_evidence_ncfg.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
   
    FL = open(path_obj["logs"] + "/glycan_citations_ncfg.log", "w")

    newrow = ["glytoucan_ac","title","journal_name","publication_date", "authors",
            "xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")]
        xref_id = row[f_list.index("evidence")]
        xref_key = "glycan_xref_pubmed"
        src_xref_id = "GLY_000528"
        src_xref_key = "glycan_xref_glygen_ds"
        combo_id = "%s %s" % (glytoucan_ac, xref_id)
        if glytoucan_ac not in glycan_list:
            FL.write("%s\n" %(glytoucan_ac))
        citerow = libgly.get_citation(xref_id, path_obj["downloads"] + "ncbi/medline/")
        if citerow != [] and combo_id not in seen:
            newrow = [glytoucan_ac] + citerow + [xref_key,xref_id, src_xref_key,src_xref_id]
            print "\"%s\"" % ("\",\"".join(newrow))
            seen[combo_id] = True

    FL.close()


    return




def extract_citations_glytoucan_ds():

    black_list = get_blacklisted_pmids()
    glycan_list = load_glycan_masterlist()

    in_file = path_obj["downloads"] + "glytoucan/current/export/pubs.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    newrow = ["glytoucan_ac","title","journal_name","publication_date", "authors",
            "xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        xref_key = "protein_xref_pubmed"
        xref_id = row[f_list.index("PubMedID")]
        source = row[f_list.index("Source")]
        source_id = row[f_list.index("SourceID")]
        src_xref_key = "glycan_xref_glytoucan"
        src_xref_id = glytoucan_ac
        if source.lower() == "unicarbkb":
            src_xref_id = source_id
            src_xref_key = "glycan_xref_unicarbkb" 
            if src_xref_id.lower().find("comp_") != -1:
                src_xref_key = "glycan_xref_unicarbkb_comp"

        original_source_id = ""
        if source_id.find("comp_") != -1:
            original_source_id = source_id
            source_id = ""
        if xref_id in black_list:
            continue
        
        combo_id = "%s %s" % (glytoucan_ac, xref_id) 
        cond_list = [glytoucan_ac not in glycan_list and glytoucan_ac != ""]
        cond_list.append(xref_id in ["0"])
        cond_list.append(combo_id in seen)
        if True in cond_list:
            continue
        citerow = libgly.get_citation(xref_id, path_obj["downloads"] + "ncbi/medline/")
        if citerow != []:
            newrow = [glytoucan_ac] + citerow + [xref_key,xref_id, src_xref_key,src_xref_id]
            print "\"%s\"" % ("\",\"".join(newrow))
            seen[combo_id] = True


    return


def extract_glytoucanidlist_ds():


    newrow = ["glytoucan_ac"]
    print "\"%s\"" % ("\",\"".join(newrow))
    data_frame = {}
    in_file = path_obj["downloads"] + "glycan_list/current/glytoucan_allacc.txt"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    newrow = [f_list[0]]
    print "\"%s\"" % ("\",\"".join(newrow))
    for row in data_frame["data"]:
        newrow = [row[0]]
        print "\"%s\"" % ("\",\"".join(newrow))

    return


def extract_masterlist_ds():

    newrow = [
        "glytoucan_ac","glytoucan_type","glycan_mass", "glycan_permass",
        "base_composition","composition","topology","monosaccharides", "is_motif",
        "missing_score"
    ]
    print "\"%s\"" % ("\",\"".join(newrow))


    missing_score_dict = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/subsumption.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        relationship = row[f_list.index("Relationship")]
        related_accession = row[f_list.index("RelatedAccession")]
        if relationship == "MissingScore":
            missing_score_dict[glytoucan_ac] = related_accession


    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/glycan_properties.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    prop_dict = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_acc")]
        newrow = [
            row[f_list.index("glytoucan_acc")],
            row[f_list.index("glytoucan_type")],
            row[f_list.index("glycan_mass")],
            row[f_list.index("glycan_permass")],
            row[f_list.index("base_composition")],
            row[f_list.index("composition")],
            row[f_list.index("topology")],
            row[f_list.index("monosaccharides")]
        ]
        prop_dict[ac] = newrow



    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/allmotifs.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    is_motif = {}
    for row in data_frame["data"]:
        motif_ac = row[f_list.index("MotifAccession")]
        prop = row[f_list.index("Property")]
        value = row[f_list.index("Value")]
        if prop == "MotifGlyTouCan":
            is_motif[value] = True
    

    img_dir = path_obj["downloads"] + "glytoucan/current/export/snfg/extended/png/"
    FL = open("logs/glycan_missing_images.log", "w")
    FL.write("glytoucan_ac,is_motif\n")

    seen_list = []
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/allglycan.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        if ac not in seen_list:
            seen_list.append(ac)
            missing_score = missing_score_dict[ac] if ac in missing_score_dict else ""
            emptyrow = [ac, "", "", "", "", "", "", "", ""]
            newrow = prop_dict[ac] if ac in prop_dict else emptyrow
            newrow = newrow + [str(ac in is_motif), missing_score] 
            print "\"%s\"" % ("\",\"".join(newrow))
            #Check image files
            img_file = img_dir + "%s.png" % (ac)
            if os.path.isfile(img_file) == False:
                FL.write("%s,%s\n" % (ac, ac in is_motif))

    FL.close()


    return



def extract_fully_determined_ds():

    glycan_list = load_glycan_masterlist()
    newrow = ["glytoucan_ac"]
    print "\"%s\""  % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/subsumption.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    seen = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        if row[f_list.index("Relationship")].strip() == "MissingScore":
            if ac not in seen and row[f_list.index("RelatedAccession")].strip() == "0":
                newrow = [ac]
                print "\"%s\""  % ("\",\"".join(newrow))
                seen[ac] = True



    return

def extract_type_n_linked_byonic_ds():

    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","byonic"]
    print "\"%s\""  % ("\",\"".join(newrow))
        
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/byonic_glygen_human_nlinked.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        bionic = row[f_list.index("Byonic")]
        if ac in glycan_list:
            newrow = [ac,bionic]
            print "\"%s\""  % ("\",\"".join(newrow))

    return


def extract_pathway_glycotree_ds():

    glycan_list = load_glycan_masterlist()
    newrow = ["glytoucan_ac","source", "target", "residue_id", "residue_name", "enzyme_uniprotkb_ac", "enzyme_tax_name"]


    print "\"%s\""  % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "/sandbox/current/allPaths.json"
    doc = json.loads(open(in_file, "r").read())
    
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")


    for obj in doc["pathways"]:
        ac = obj["glytoucan_ac"]
        if ac in glycan_list and "links" in obj:
            for o in obj["links"]:
                source, target = o["source"], o["target"]
                if o["residue_affected"] == []:
                    continue
                residue_id = o["residue_affected"]["residue_id"]
                residue_name = o["residue_affected"]["full_name"]
                enzyme_list = []
                if "enzymes" in o:
                    for oo in o["enzymes"]:
                        species = oo["species"]
                        tax_name = species_obj[species]["long_name"] if species in species_obj else species
                        enzyme_list.append({"ac":oo["uniprot"], "tax_name":tax_name})
                enzyme_list = [{"ac":"", "tax_name":""}] if enzyme_list == [] else enzyme_list
                for oo in enzyme_list:
                    newrow = [ac, source, target, residue_id, residue_name, oo["ac"], oo["tax_name"]]
                    print "\"%s\""  % ("\",\"".join(newrow))


    return


def extract_dictionary_ds():


    #term (main_entry),glycan_dictionary_accession,glytoucan_accession ,
    #term_in_sentence,publication,definition,term_xref,best_match,synonymns,
    #function,disease_associations,wikipedia,essentials_of_glycobiology

    glycan_list = load_glycan_masterlist()


    data_frame = {}
    in_file = path_obj["downloads"] + "glycan_dictionary/current/glycan_dictionary.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
   
    ignore_fields = [
        "term (main_entry)", 
        "glycan_dictionary_accession", 
        "glytoucan_accession"
    ]
    #ignore_fields = []

    newrow = ["glytoucan_ac","term"]
    for f in f_list:
        if f in ignore_fields:
            continue
        newrow.append(f)
    newrow += ["xref_key","xref_id"]

    print "\"%s\""  % ("\",\"".join(newrow))

    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_accession")]
        #if ac not in glycan_list:
        #    continue
        term = row[0]
        glycan_dictionary_accession = row[f_list.index("glycan_dictionary_accession")]
        newrow = [ac, term]
        for f in f_list:
            if f in ignore_fields:
                continue
            newrow.append(row[f_list.index(f)].replace("\"", "\'"))
        newrow += ["glycan_xref_dictionary", glycan_dictionary_accession]
        print "\"%s\""  % ("\",\"".join(newrow))


    return 

def extract_species_ds():

    glycan_list = load_glycan_masterlist()
    motif_list = load_motif_masterlist()

    newrow = ["glytoucan_ac","tax_id","tax_name","annotation_category","source","source_id",
            "original_source_id", "xref_key", "xref_id","is_motif"]
    print "\"%s\""  % ("\",\"".join(newrow))


    seen = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/species_expanded.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    
    sourceid2taxid = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCan AC")]
        source_id = row[f_list.index("Source ID")]
        tax_id = row[f_list.index("tax_id")]
        cat = row[f_list.index("AnnotationCategory")]
        if cat == "Direct":
            #if cat == Direct then source_id = glytoucan_ac
            source_id = glytoucan_ac
            if source_id not in sourceid2taxid:
                sourceid2taxid[source_id] = []
            if tax_id.strip() != "" and tax_id not in sourceid2taxid[source_id]:
                sourceid2taxid[source_id].append(tax_id)

    seen = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCan AC")]
        species = row[f_list.index("Species Name")]
        cat = row[f_list.index("AnnotationCategory")]
        source = row[f_list.index("Source")]
        source_id = row[f_list.index("Source ID")]
        tax_id = row[f_list.index("tax_id")]
        if tax_id == "11108":
            tax_id = "63746"

        original_source_id = ""
        species = "sarscov2" if species == "SARS" else species
        is_motif = "True" if ac in motif_list else "False"
        if source_id.lower().find("comp_") != -1:
            original_source_id = source_id
            source_id = "UniCarbKB"
        if cat == "Direct":
            tax_name = species_obj[tax_id]["long_name"]
            xref_id = source_id
            xref_key = "glycan_xref_" + source.lower()
            xref_key = "glycan_xref_glygen_ds" if xref_key == "glycan_xref_glygen" else xref_key
            xref_key = "glycan_xref_unicarbkb_comp" if source_id == "UniCarbKB" else xref_key
            newrow = [ac,tax_id,tax_name, cat, source, source_id,original_source_id,
                xref_key,xref_id]
            row_str = json.dumps(newrow)
            if row_str not in seen:
                is_ref = species_obj[tax_id]["is_reference"] if tax_id in species_obj else "no"
                #print tax_id, is_ref
                print "\"%s\""  % ("\",\"".join(newrow + [is_motif]))
                seen[row_str] = True
        elif cat  == "Subsumption":
            xref_id = source_id
            xref_key = "glycan_xref_" + source.lower()
            xref_key = "glycan_xref_glygen_ds" if xref_key == "glycan_xref_glygen" else xref_key
            if source_id in sourceid2taxid:
                row_list_one, row_list_two, row_list_three = [], [], []
                for tax_id in sourceid2taxid[source_id]:
                    if tax_id == "11108":
                        tax_id = "63746"
                    tax_name_one = species
                    if species.lower() in species_obj:
                        tax_name_one = species_obj[species.lower()]["long_name"]
                    tax_name_two = species_obj[tax_id]["long_name"]
                    newrow = [ac,tax_id,tax_name_two, cat, source, source_id,original_source_id,
                        xref_key,xref_id]
                    if tax_name_one == tax_name_two:
                        row_list_one.append(newrow)
                    elif tax_name_one.find(tax_name_two) != -1:
                        row_list_two.append(newrow)
                    else:
                        row_list_three.append(newrow)
                
                if row_list_one != []:          
                    for newrow in row_list_one:
                        row_str = json.dumps(newrow)
                        if row_str not in seen:
                            is_ref = species_obj[tax_id]["is_reference"] if tax_id in species_obj else "no"
                            #print tax_id, is_ref
                            print "\"%s\""  % ("\",\"".join(row_list_one[0] + [is_motif]))
                            seen[row_str] = True
                #elif row_list_two != []:
                if row_list_two != []:
                    for newrow in row_list_two:
                        row_str = json.dumps(newrow)
                        if row_str not in seen:
                            is_ref = species_obj[tax_id]["is_reference"] if tax_id in species_obj else "no"
                            #print tax_id, is_ref
                            print "\"%s\""  % ("\",\"".join(newrow + [is_motif]))
                            seen[row_str] = True
                #else:
                if row_list_three != []:
                    for newrow in row_list_three:
                        row_str = json.dumps(newrow)
                        if row_str not in seen:
                            is_ref = species_obj[tax_id]["is_reference"] if tax_id in species_obj else "no"
                            #print tax_id, is_ref
                            print "\"%s\""  % ("\",\"".join(newrow + [is_motif]))
                            seen[row_str] = True

    return



def extract_species_ds_old():


    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","tax_id","tax_name","annotation_category","source","source_id",
            "original_source_id", "xref_key", "xref_id"]
    print "\"%s\""  % ("\",\"".join(newrow))


    seen = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/species.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    #load ann_cat_dict and ann_evd_dict
    ann_cat_dict = {}
    ann_evd_dict = {}


    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        species = row[f_list.index("Species")]
        value = row[f_list.index("Value")]
        if species.find("Category") != -1:
            sp_name, cat_value = species.split(" ")[0].lower(), value
            sp_name = "sarscov2" if sp_name == "sars" else sp_name
            if ac not in ann_cat_dict:
                ann_cat_dict[ac] = {}
            if sp_name not in ann_cat_dict[ac]:
                ann_cat_dict[ac][sp_name] = []
            ann_cat_dict[ac][sp_name].append(cat_value)
        elif species.find("Evidence") != -1:
            sp_name, evd_value = species.split(" ")[0].lower(), value
            sp_name = "sarscov2" if sp_name == "sars" else sp_name
            if ac not in ann_evd_dict:
                ann_evd_dict[ac] = {}
            if sp_name not in ann_evd_dict[ac]:
                ann_evd_dict[ac][sp_name] = []
            ann_evd_dict[ac][sp_name].append(evd_value)


    out_list = []
    ac_list = list(set(ann_evd_dict.keys() + ann_cat_dict.keys()))
    for ac in ac_list:
        sp_list = ann_cat_dict[ac].keys() if ac in ann_cat_dict else []
        sp_list += ann_evd_dict[ac].keys() if ac in ann_evd_dict else []
        for sp_name in list(set(sp_list)):
            cat_list = []
            if ac in ann_cat_dict:
                cat_list = ann_cat_dict[ac][sp_name] if sp_name in ann_cat_dict[ac] else []
            
            evd_list = []
            if ac in ann_evd_dict:
                evd_list = ann_evd_dict[ac][sp_name] if sp_name in ann_evd_dict[ac] else []

            
            for cat in cat_list:
                if cat == "Subsumption" and sp_name not in ["rat", "hcv"]:
                    for evd in evd_list:
                        #tax_id = ""
                        tax_id = str(species_obj[sp_name]["tax_id"])
                        tax_name = species_obj[sp_name]["long_name"]
                        source = "Subsumption"
                        xref_id = evd.split(" ")[-1]
                        xref_key = "glycan_xref_subsumption"
                        source_id = "via " + evd.split(" ")[-1]
                        newrow = [ac,tax_id,tax_name,cat,source,source_id,"",xref_key,xref_id]
                        if ac in glycan_list:
                            newrow_str = ",".join(newrow)
                            if newrow_str not in seen:
                                out_list.append(newrow)
                                seen[newrow_str] = True
                elif cat == "Composition":
                    for evd in evd_list:
                        tax_id = ""
                        tax_name = species_obj[sp_name]["long_name"]
                        source = "Composition"
                        xref_id = evd.split(" ")[-1]
                        xref_key = "glycan_xref_composition"
                        source_id = "via " + evd.split(" ")[-1]
                        newrow = [ac,tax_id,tax_name,cat,source,source_id,"",xref_key,xref_id]
                        if ac in glycan_list:
                            newrow_str = ",".join(newrow)
                            if newrow_str not in seen:
                                #Stop outputing these rows (Rahi)
                                #print "\"%s\""  % ("\",\"".join(newrow))
                                seen[newrow_str] = True
                elif cat == "Direct":
                    for evd in evd_list:
                        if evd.find("GlyTouCan") != -1:
                            source = "GlyTouCan"
                            source_id = ac
                            xref_id = ac
                            xref_key = "glycan_xref_glytoucan"
                            tax_id = evd.split(" ")[-1].strip()
                            tax_name = species_obj[tax_id]["long_name"]
                            newrow = [ac,tax_id, tax_name,cat, source, source_id, "",
                                    xref_key, xref_id]
                            if ac in glycan_list:
                                newrow_str = ",".join(newrow)
                                if newrow_str not in seen:
                                    out_list.append(newrow) 
                                    seen[newrow_str] = True
                        elif evd.find("UniCarbKB") != -1:
                            source = "UniCarbKB"
                            source_id = evd.split(" ")[-3].split(":")[1]
                            xref_id = source_id
                            xref_key = "glycan_xref_unicarbkb"
                            tax_id = evd.split(" ")[-1].strip()
                            tax_name = ""
                            if tax_id in species_obj:
                                tax_name = species_obj[tax_id]["long_name"]
                            if source_id.lower().find("comp_") == -1:
                                original_source_id = ""
                                newrow = [ac,tax_id, tax_name, cat, source, source_id, 
                                        original_source_id, xref_key,xref_id]
                                if ac in glycan_list:
                                    newrow_str = ",".join(newrow)
                                    if newrow_str not in seen:
                                        out_list.append(newrow) 
                                        seen[newrow_str] = True
                            else:
                                original_source_id = source_id
                                source_id = ""
                                xref_key = "glycan_xref_unicarbkb_comp"
                                xref_id = source_id
                                tax_name = ""
                                if tax_id in species_obj:
                                    tax_name = species_obj[tax_id]["long_name"]
                                newrow = [ac,tax_id, tax_name, cat, source, source_id, 
                                        original_source_id, xref_key,xref_id]
                                if ac in glycan_list:
                                    newrow_str = ",".join(newrow)
                                    if newrow_str not in seen:
                                        out_list.append(newrow) 
                                        seen[newrow_str] = True
                        #elif evd.find(":") != -1:
                        elif evd.find("TaxID") != -1:
                            source = evd.split(" ")[-3].split(":")[0]
                            source_id = ""
                            if evd.find(":") != -1:
                                source_id = evd.split(" ")[-3].split(":")[1]
                            #original_source_id = source_id
                            original_source_id = ""
                            xref_key = "glycan_xref_" + source.lower()
                            xref_id = source_id
                            tax_id = evd.split(" ")[-1].strip()
                            tax_name = ""
                            if tax_id in species_obj:
                                tax_name = species_obj[tax_id]["long_name"]
                            newrow = [ac,tax_id, tax_name, cat, source, source_id,
                                        original_source_id, xref_key,xref_id]
                            if ac in glycan_list:
                                newrow_str = ",".join(newrow)
                                if newrow_str not in seen:
                                    out_list.append(newrow) 
                                    seen[newrow_str] = True

    for newrow in out_list:
        print "\"%s\""  % ("\",\"".join(newrow))
    return




    direct_ann_dict = {}
    for row in out_list:
        ac, tax_id, tax_name, ann_type = row[0], row[1], row[2], row[3]
        if ann_type == "Direct":
            if ac not in direct_ann_dict:
                direct_ann_dict[ac] = {}
            direct_ann_dict[ac][tax_id] = tax_name
    
    seen_row = {}
    for row in out_list:
        ac, tax_id, tax_name, ann_type = row[0], row[1], row[2], row[3]
        if ann_type == "Subsumption":
            s_ac = row[5].split(" ")[1]
            for s_tax_id in direct_ann_dict[s_ac]:
                newrow = []
                for v in row:
                    newrow.append(v)
                newrow[1] = s_tax_id
                newrow[2] = direct_ann_dict[s_ac][s_tax_id]
                row_str = json.dumps(newrow)
                if row_str not in seen_row:
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen_row[row_str] = True
        else:
            newrow = row
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                print "\"%s\""  % ("\",\"".join(newrow))
                seen_row[row_str] = True



    return


def extract_compiled_ds(ds_name, filter_flag):

    glycan_list = []
    if filter_flag == True:
        glycan_list = load_glycan_masterlist()

    in_file = "compiled/glycan_%s.csv" % (ds_name)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    newrow = f_list
    if ds_name == "synthesized":
        #newrow[f_list.index("evidence")] = "xref_id"
        newrow += ["xref_key", "xref_id"]
    print "\"%s\""  % ("\",\"".join(newrow))
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        row_list = []
        if ds_name == "synthesized":
            xref_id = row[f_list.index("evidence")]
            xref_key = "protein_xref_pubmed"
            row_list.append(row + [xref_key, xref_id])
            row_list.append(row + ["glycan_xref_glygen_ds", "GLY_000309"])
        else:
            row_list.append(row)

        if filter_flag == True:
            if ac in glycan_list:
                for newrow in row_list:
                    print "\"%s\""  % ("\",\"".join(newrow))
        else:
            for newrow in row_list:
                print "\"%s\""  % ("\",\"".join(newrow))

    return

def extract_evidence_ncfg_ds():

    glycan_list = load_glycan_masterlist()

    in_file = path_obj["downloads"] + "ncfg/current/glycan_evidence_ncfg.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    newrow = f_list
    print "\"%s\""  % ("\",\"".join(newrow))
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        newrow = row
        if ac in glycan_list:
            print "\"%s\""  % ("\",\"".join(newrow))

    return


def extract_xrefs_ds(dataset):

    filename_dict = {
        "xref_pubchem":"pbch"
    }

    ds_filename = "_".join(dataset.split("_")[1:])
    ds_suffix = "glycan_" + dataset
    if dataset in filename_dict:
        ds_filename = filename_dict[dataset]

    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","xref_id","xref_key"]
    print "\"%s\""  % ("\",\"".join(newrow))

    seen_row = {} 
    if dataset == "xref_glytoucan":
        data_frame = {}
        in_file = path_obj["unreviewed"] + "glycan_masterlist.csv"
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            glytoucan_ac = row[f_list.index("glytoucan_ac")]
            newrow = [glytoucan_ac, glytoucan_ac, "glycan_xref_glytoucan"]
            print "\"%s\""  % ("\",\"".join(newrow))
    elif dataset == "xref_rhea":
        cid2rheaid = {}
        in_file = path_obj["downloads"] + "chebi/current/reference.tsv"
        cmd = "head -1 " + in_file
        f_list = commands.getoutput(cmd).strip().split("\t")
        cmd = "grep RHEA: " + in_file
        lines = commands.getoutput(cmd).split("\n")
        for line in lines:
            row = line.split("\t")
            c_id = row[f_list.index("COMPOUND_ID")]
            rhea_id = row[f_list.index("REFERENCE_ID")].split(":")[1]
            if c_id not in cid2rheaid:
                cid2rheaid[c_id] = {}
            cid2rheaid[c_id][rhea_id] = True
        data_frame = {}
        in_file = path_obj["unreviewed"] + "glycan_xref_chebi.csv"
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            glytoucan_ac = row[f_list.index("glytoucan_ac")]
            c_id = row[f_list.index("xref_id")]
            if c_id in cid2rheaid:
                for rhea_id in cid2rheaid[c_id]:
                    newrow = [glytoucan_ac, rhea_id, "glycan_xref_rhea"]
                    row_str = json.dumps(newrow)
                    if row_str not in seen_row:
                        print "\"%s\""  % ("\",\"".join(newrow))
                        seen_row[row_str] = True

    elif dataset == "xref_gptwiki":
        data_frame = {}
        in_file = path_obj["downloads"] + "/gptwiki/current/glycosites.csv"
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            glytoucan_ac = row[f_list.index("GlyTouCan")]
            if glytoucan_ac in glycan_list:
                newrow = [glytoucan_ac, glytoucan_ac, "glycan_xref_gptwiki"]
                row_str = json.dumps(newrow)
                if row_str not in seen_row:
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen_row[row_str] = True
    elif dataset == "xref_dictionary":
        data_frame = {}
        in_file = path_obj["downloads"] + "/glycan_dictionary/current/glycan_dictionary.csv"
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            glytoucan_ac = row[f_list.index("glytoucan_accession")]
            glycan_dictionary_accession = row[f_list.index("glycan_dictionary_accession")]
            if glytoucan_ac in glycan_list:
                newrow = [glytoucan_ac,glycan_dictionary_accession, "glycan_xref_dictionary"]
                row_str = json.dumps(newrow)
                if row_str not in seen_row:
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen_row[row_str] = True
    elif dataset == "xref_reactome":
        cid2reactomeid = {}
        data_frame = {}
        in_file = path_obj["downloads"] + "reactome/current/ChEBI2Reactome_PE_Pathway.csv"
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            c_id = row[f_list.index("chebi_id")]
            reactome_id = row[f_list.index("reactome_compound_id")]
            if c_id not in cid2reactomeid:
                cid2reactomeid[c_id] = {}
            cid2reactomeid[c_id][reactome_id] = True
        data_frame = {}
        in_file = path_obj["unreviewed"] + "glycan_xref_chebi.csv"
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            glytoucan_ac = row[f_list.index("glytoucan_ac")]
            c_id = row[f_list.index("xref_id")]
            if c_id in cid2reactomeid:
                for reactome_id in cid2reactomeid[c_id]:
                    newrow = [glytoucan_ac, reactome_id, "glycan_xref_reactome"]
                    row_str = json.dumps(newrow)
                    if row_str not in seen_row:
                        print "\"%s\""  % ("\",\"".join(newrow))
                        seen_row[row_str] = True

    else:
        data_frame = {}
        if dataset == "xref_sandbox":
            #in_file = path_obj["unreviewed"] + "/glycan_enzyme.csv"
            in_file = path_obj["downloads"] + "sandbox/current/sandbox_accessions.csv"
            libgly.load_sheet(data_frame, in_file, ",")
        else:
            in_file = path_obj["downloads"] + "glytoucan/current/export/%s.tsv" % (ds_filename)
            libgly.load_sheet(data_frame, in_file, "\t")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[0]
            xref_id = ac if dataset in ["xref_sandbox", "xref_gnome", "xref_glycosmos"] else row[1]
            if xref_id == "":
                continue
            if ac not in glycan_list:
                continue
            xref_key = "glycan_" + dataset
            if dataset == "xref_gnome":
                xref_key = "glycan_xref_gnome_cross_reference"
            if dataset == "xref_pubchem":
                xref_key = "glycan_xref_pubchem_compound" if xref_id[0:3] == "CID" else "glycan_xref_pubchem_substance"
            if dataset == "xref_unicarbkb":
                if xref_id.lower().find("comp_") != -1 or xref_id.lower().find("hex") != -1:
                    xref_key = "glycan_xref_unicarbkb_comp"
            
            if dataset == "xref_glyconnect":
                ac_type = row[f_list.index("GlyConnectAccessionType")].lower()
                xref_key = "glycan_xref_" + ac_type
            
            xref_id = xref_id.replace("CID", "")
            xref_id = xref_id.replace("SID", "")
            newrow = [ac, xref_id, xref_key]
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                #print "xxxxx"
                print "\"%s\""  % ("\",\"".join(newrow))
                seen_row[row_str] = True
                
    return



def extract_enzyme_ds():

    glycan_list = load_glycan_masterlist()

    canon2genename = {}
    ac2canon = {}
    file_list = glob.glob(path_obj["unreviewed"] +  "*_protein_masterlist.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            canon2genename[canon] = row[f_list.index("gene_name")]
            for isoform in [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]:
                ac = isoform.split("-")[0]
                ac2canon[ac] = canon

    canon2recname = {}
    file_list = glob.glob(path_obj["unreviewed"] +  "*_protein_recnames.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            canon2recname[canon] = row[f_list.index("recommended_name_full")]

    newrow = ["glytoucan_ac","residue_name","residue_id","uniprotkb_canonical_ac","gene_name","gene_id","parent_residue_id","enzyme_type","species","recommended_name_full", "xref_key", "xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))

    data_frame = {}
    
    in_file = path_obj["downloads"] + "sandbox/current/glycotree_annotated_glycans.tsv"
    #libgly.load_sheet(data_frame, in_file, ",")
    libgly.load_sheet(data_frame, in_file, "\t")


    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        val_dict = {}
        for f in f_list:
            val_dict[f] = row[f_list.index(f)]
        if val_dict["glytoucan_ac"] not in glycan_list:
            continue
        if val_dict["uniprot"] not in ac2canon:
            continue
        canon = ac2canon[val_dict["uniprot"]]
        recname = canon2recname[canon] if canon in canon2recname else ""
        newrow = [
            val_dict["glytoucan_ac"],val_dict["residue_name"],val_dict["residue_id"],canon,
            val_dict["gene_name"],val_dict["gene_id"],val_dict["parent_residue_id"],val_dict["enzyme_type"],
            val_dict["species"], recname,"glycan_xref_glygen_ds", "GLY_000284"
        ]
        print "\"%s\"" % ("\",\"".join(newrow))

    return





def get_glyco_stats(species_list):

    used_list = get_used_ds_list()
    file_list = glob.glob(path_obj["unreviewed"] + "/*glycosylation_sites_*.csv")
    seen = {"glycan":{}, "glycanpos":{}, "pos":{}, "siteinfo":{}}
    for in_file in file_list:
        file_name = in_file.split("/")[-1]
        if file_name not in used_list:
            continue
        species = in_file.split("/")[-1].split("_")[0]
        if species_list != [] and species not in species_list:
            continue
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            glycan = row[f_list.index("saccharide")]
            amino_acid = row[f_list.index("amino_acid")]
            gly_type = row[f_list.index("glycosylation_type")].lower()
            aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]

            if canon not in seen["pos"]:
                seen["pos"][canon] = {}
            if aa_pos not in seen["pos"][canon]:
                seen["pos"][canon][aa_pos] = True
           
            if xref_key == "protein_xref_pubmed":
                site_combo = "%s|%s" % (amino_acid,gly_type)
                if canon not in seen["siteinfo"]:
                    seen["siteinfo"][canon] = {}
                if aa_pos not in seen["siteinfo"][canon]:
                    seen["siteinfo"][canon][aa_pos] = {}
                if site_combo not in seen["siteinfo"][canon][aa_pos]:
                    seen["siteinfo"][canon][aa_pos][site_combo] = {}
                seen["siteinfo"][canon][aa_pos][site_combo][xref_id] = True
            

            #if canon == "P00533-1":
            #    print "flag", aa_pos, amino_acid, gly_type, glycan

            if canon not in seen["glycan"]:
                seen["glycan"][canon] = {}
            if gly_type not in seen["glycan"][canon]:
                seen["glycan"][canon][gly_type] = {}
            if glycan != "" and glycan not in seen["glycan"][canon][gly_type]:
                seen["glycan"][canon][gly_type][glycan] = True
            
            if canon not in seen["glycanpos"]:
                seen["glycanpos"][canon] = {}
            if gly_type not in seen["glycanpos"][canon]:
                seen["glycanpos"][canon][gly_type] = {}
            if glycan != "" and glycan not in seen["glycanpos"][canon][gly_type]:
                seen["glycanpos"][canon][gly_type][aa_pos] = True


    return seen


def get_used_ds_list():
    
    file_name_list = []
    ds_obj_list = json.loads(open("generated/misc/dataset-masterlist.json", "r").read())
    for obj in ds_obj_list:
        ds_name = obj["name"]
        ds_format = obj["format"]
        mol = obj["categories"]["molecule"]
        if ds_format not in ["csv"]:
            continue
        if obj["categories"]["species"] == []:
            if obj["integration_status"]["status"] == "integrate_all":
                if "protein" in obj["target_objects"]:
                    file_name_list.append("%s_%s.csv" % (mol, ds_name))
        elif obj["integration_status"]["status"] != "integrate_none":
            sp_list_one = sorted(obj["categories"]["species"])
            for species in sp_list_one:
                if species not in obj["integration_status"]["excludelist"]:
                    if "protein" in obj["target_objects"]:
                        file_name_list.append("%s_%s_%s.csv" % (species, mol, ds_name))


    return file_name_list



def get_used_ds_list_old():     

    species_list = []
    for k in species_obj:
        obj = species_obj[k]
        if obj["short_name"] not in species_list and obj["is_reference"] == "yes":
            species_list.append(obj["short_name"])


    seen_ds = {}
    in_file = "generated/misc/protein_datasets.json"
    file_list_obj = json.loads(open(in_file, "r").read())
    for cat in file_list_obj:
        if cat == "common":
            for species in species_list:
                for mol in file_list_obj[cat]:
                    for sheet_name in file_list_obj[cat][mol]:
                        ds = "%s_%s_%s.csv" % (species, mol, sheet_name)
                        seen_ds[ds] = True
        else:
            for mol in file_list_obj[cat]:
                for sheet_name in file_list_obj[cat][mol]:
                    ds = "%s_%s_%s.csv" % (cat, mol, sheet_name)
                    if cat == "agnostic":
                        ds = "%s_%s.csv" % (mol, sheet_name)
                    seen_ds[ds] = True

    in_file = "generated/misc/glycan_datasets.json"
    file_list_obj = json.loads(open(in_file, "r").read())
    for sheet_name in file_list_obj["glycan"]:
        ds = "glycan_%s.csv" % (sheet_name)
        seen_ds[ds] = True

    return seen_ds.keys()





def extract_xref_mapping_ds(ds_name):
   
    map_dict = {
        "glygen_genecards_xref_mapping":{
            "headers":["uniprotkb_primary_accession","additional_info","glygen_url"]
        },
        "glygen_uniprotkb_xref_mapping":{
            "headers":["uniprotkb_primary_accession","additional_info"]
        },
        "glygen_pharos_xref_mapping":{
            "headers":["uniprotkb_ac","glygen_protein_url"]
        },
        "glygen_pubchem_xref_mapping":{
            "headers":["uniprotkb_ac","refseq_ac", "glycosylation_annotation"]
        },
        "glygen_iptmnet_xref_mapping":{
            "headers":["uniprotkb_ac", "glycosylation_site_uniprotkb", "amino_acid",
                "glycosylation_type", "pmid", "glygen_url"
            ]
        }
    }

    gtc2cid = {}    
    in_file = path_obj["unreviewed"] + "/glycan_xref_pubchem.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        gtc = row[f_list.index("glytoucan_ac")]
        xref_id = row[f_list.index("xref_id")]
        xref_key = row[f_list.index("xref_key")]
        if xref_key == "glycan_xref_pubchem_compound":
            if gtc not in gtc2cid:
                gtc2cid[gtc] = []
            gtc2cid[gtc].append(xref_id)

    glycan_site_dict = {}
    file_list = glob.glob(path_obj["unreviewed"] + "/*_proteoform_glycosylation_sites*.csv")
    for in_file in file_list:
        species = in_file.split("/")[-1].split("_")[0]
        tax_id = species_obj[species]["tax_id"]
        tax_name = species_obj[species]["long_name"]
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            saccharide = row[f_list.index("saccharide")]
            xref_key = row[f_list.index("xref_key")]
            if saccharide == "":
                continue
            if xref_key not in ["protein_xref_pubmed"]:
                continue
            amino_acid = row[f_list.index("amino_acid")]
            aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
            gly_type = row[f_list.index("glycosylation_type")]
            pmid = row[f_list.index("xref_id")]
            if canon not in glycan_site_dict:
                glycan_site_dict[canon] = []
            if canon not in glycan_site_dict:
                glycan_site_dict[canon] = []
            o = {"amino_acid":amino_acid, "aa_pos":aa_pos, "saccharide":saccharide, 
                    "gly_type":gly_type, "pmid":pmid}
            glycan_site_dict[canon].append(o)



    uniprotkbac2refseqac = {}
    file_list = glob.glob(path_obj["unreviewed"] + "/*_protein_xref_refseq.csv")
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            uniprotkb_ac = row[f_list.index("uniprotkb_canonical_ac")].split("-")[0]
            refseq_ac = row[f_list.index("xref_id")]
            uniprotkbac2refseqac[uniprotkb_ac] = refseq_ac


    canon_dict = {}
    ac2canon = {}
    file_list = glob.glob(path_obj["unreviewed"] + "/*_protein_masterlist.csv")
    for in_file in file_list:
        species = in_file.split("/")[-1].split("_")[0]
        tax_id = species_obj[species]["tax_id"]
        tax_name = species_obj[species]["long_name"]
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            ac = canon.split("-")[0]
            tax_name = species_obj[species]["long_name"]
            canon_dict[canon] = {"ac":ac, "taxid":tax_id, "taxname":tax_name}
            ac2canon[ac] = canon

    seen_row = {}
    if ds_name in ["glygen_uniprotkb_xref_mapping", "glygen_genecards_xref_mapping", "glygen_pubchem_xref_mapping"]:
        sp_list = ["human"] if ds_name == "glygen_genecards_xref_mapping" else []
        glyco_stats = get_glyco_stats(sp_list)

        newrow = map_dict[ds_name]["headers"]
        if ds_name == "glygen_uniprotkb_xref_mapping":
            print "%s" % ("\t".join(newrow))
        elif ds_name == "glygen_genecards_xref_mapping": 
            print "%s"  % ("\t".join(newrow))
        elif ds_name == "glygen_pubchem_xref_mapping":
            newrow += ["amino_acid","glycosylation_site_uniprotkb",
                    "saccharide","pubchem_cid","glycosylation_type","pmid"]
            print "\"%s\"" % ("\",\"".join(newrow))
        for canon in glyco_stats["glycan"]:
            row = []
            total_n = 0
            for gly_type in sorted(glyco_stats["glycan"][canon]):
                n = len(glyco_stats["glycan"][canon][gly_type].keys())
                n_g = len(glyco_stats["glycanpos"][canon][gly_type].keys())
                if n > 0:
                    sf_one = "glycans" if n > 1 else "glycan"
                    sf_two = "sites" if n_g > 1 else "site"
                    row.append("%s %s %s (%s %s)" % (n, gly_type[0].upper() + gly_type[1:], 
                        sf_one, n_g, sf_two))
                    total_n += n
            n_pos = len(glyco_stats["pos"][canon].keys()) if canon in glyco_stats["pos"] else 0
            ac = canon.split("-")[0]
            suffix = "sites" if n_pos > 1 else "site"
            newrow = [ac, "%s %s, %s" % (n_pos, suffix, ", ".join(row))]

            if total_n == 0:
                suffix = "sites" if n_pos > 1 else "site"
                newrow = [ac,  "%s %s" % (n_pos, suffix)]
            if ds_name == "glygen_pubchem_xref_mapping":
                refseq_ac = uniprotkbac2refseqac[newrow[0]] if newrow[0] in uniprotkbac2refseqac else ""
                newrow = [newrow[0], refseq_ac, newrow[1]]
                if canon not in glycan_site_dict:
                    newnewrow = newrow +  ["","","","","",""]
                    row_str = json.dumps(newnewrow)
                    if row_str not in seen_row:
                        print "\"%s\"" % ("\",\"".join(newnewrow))
                        seen_row[row_str] = True
                else:
                    
                    for o in glycan_site_dict[canon]:
                        gtc = o["saccharide"]
                        if gtc not in gtc2cid:
                            continue
                        for cid in gtc2cid[gtc]:
                            newnewrow = newrow +  [o["amino_acid"],o["aa_pos"],o["saccharide"],
                                    cid,o["gly_type"],o["pmid"]]
                            row_str = json.dumps(newnewrow)
                            if row_str not in seen_row:
                                #print "Robel-3", canon, gtc, cid
                                print "\"%s\"" % ("\",\"".join(newnewrow))
                                seen_row[row_str] = True
            elif ds_name == "glygen_genecards_xref_mapping":
                glygen_url = "https://glygen.org/protein/%s" % (newrow[0])
                newrow.append(glygen_url)
                print "%s" % ("\t".join(newrow))
            else:
                print "%s" % ("\t".join(newrow))
    elif ds_name == "glygen_iptmnet_xref_mapping":
        glyco_stats = get_glyco_stats(["human", "mouse", "rat"])
        newrow = map_dict[ds_name]["headers"]
        print "\"%s\"" % ("\",\"".join(newrow))
        for canon in sorted(glyco_stats["siteinfo"]):
            ac = canon.split("-")[0]
            for pos in sorted(glyco_stats["siteinfo"][canon]):
                if pos == "":
                    continue
                glygen_url = "https://www.glygen.org/Siteview/%s/%s" % (ac, pos)
                for site_combo in glyco_stats["siteinfo"][canon][pos]:
                    amino_acid,gly_type = site_combo.split("|")
                    pmid = "|".join(glyco_stats["siteinfo"][canon][pos][site_combo].keys())
                    newrow = [ac, str(pos), amino_acid, gly_type, pmid, glygen_url]
                    print "\"%s\"" % ("\",\"".join(newrow))
    elif ds_name == "glygen_pharos_xref_mapping":
        FL = open("logs/glygen_pharos_xref_mapping.log", "w")
        FL.write("unmapped_accession\n")
        newrow = map_dict[ds_name]["headers"]
        print "\"%s\"" % ("\",\"".join(newrow))
        file_list = glob.glob(path_obj["downloads"] + "pharos/current/batch.*.json")
        for in_file in file_list:
            doc = json.loads(open(in_file, "r").read())
            for obj in doc["content"]:
                ac = obj["accession"]
                if ac in ac2canon:
                    canon = ac2canon[ac]
                    glygen_url = "https://glygen.org/protein/%s" % (ac)
                    newrow = [ac, glygen_url]
                    print "\"%s\"" % ("\",\"".join(newrow))
                else:
                    FL.write("%s\n" % (ac))
        FL.close()

    return


def extract_xref_chebi_ds():

    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","xref_id","xref_key"]
    print "\"%s\""  % ("\",\"".join(newrow))

    chebi2inchi = {}
    in_file = path_obj["downloads"] + "chebi/current/database_accession.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        source = row[f_list.index("SOURCE")]
        if source != "GlyTouCan":
            continue
        glytoucan_ac = row[f_list.index("ACCESSION_NUMBER")]
        xref_id = row[f_list.index("COMPOUND_ID")]
        xref_key = "glycan_xref_chebi"
        if glytoucan_ac not in glycan_list:
            continue
        newrow = [glytoucan_ac, xref_id, xref_key]
        print "\"%s\"" % ("\",\"".join(newrow))

    return




def extract_sequences_smiles_isomeric_ds():

    cid2glytoucan = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_xref_pubchem.csv" 
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "glycan_xref_pubchem_compound":
            cid2glytoucan[xref_id] = ac



    newrow = ["glytoucan_ac","pubchem_id", "sequence_smiles_isomeric"]
    print "\"%s\"" % ("\",\"".join(newrow))
    in_file = path_obj["downloads"] + "pubchem/compound/current/CID-SMILES"

    with open(in_file, "r") as FR:
        for line in FR:
            cid, smiles = line.strip().split("\t")
            if cid in cid2glytoucan:
                newrow = [cid2glytoucan[cid], cid, str(smiles)]
                print "\"%s\"" % ("\",\"".join(newrow))



    return




def extract_sequences_inchi_ds():

    glycan_list = load_glycan_masterlist()

    cid2glytoucan = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "glycan_xref_pubchem.csv" 
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "glycan_xref_pubchem_compound":
            cid2glytoucan[xref_id] = ac


    newrow = ["glytoucan_ac","sequence_inchi","inchi_key"]
    print "\"%s\"" % ("\",\"".join(newrow))

    data_frame = {}
    in_file = "compiled/cid2inchi.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        cid = row[f_list.index("pubchem_cid")]
        inchi = row[f_list.index("inchi")]
        inchikey = row[f_list.index("inchikey")]
        if cid in cid2glytoucan:
            glytoucan_ac = cid2glytoucan[cid]
            if glytoucan_ac in glycan_list:
                newrow = [glytoucan_ac, inchi, inchikey]
                print "\"%s\"" % ("\",\"".join(newrow))


    return


def extract_disease_names_ds():

    species_list = []
    for k in species_obj:
        if k.isdigit() == True:
            continue
        if species_obj[k]["is_reference"] == "no":
            continue
        species_list.append(k)

    
    mondo_dict = {}
    mondo2name = {}
    mondo2synonym = {}
    mondo2desc = {}
    seen_row = {}

    mondo_obj = json.loads(open(path_obj["downloads"] + "ohsu/mondo.json", "r").read())
    for g_obj in mondo_obj["graphs"]:
        for n_obj in g_obj["nodes"]:
            mondo_id = n_obj["id"].split("/")[-1]
            if mondo_id[0:6] == "MONDO_":
                m_id = mondo_id[6:]
                if "lbl" not in n_obj:
                    continue
                mondo2name[m_id] = n_obj["lbl"]
                if "meta" in n_obj:
                    if "synonyms" in n_obj["meta"]:
                        for o in n_obj["meta"]["synonyms"]:
                            if m_id not in mondo2synonym:
                                mondo2synonym[m_id] = []
                            if o["val"] not in mondo2synonym[m_id]:
                                mondo2synonym[m_id].append(o["val"])
                    if "definition" in n_obj["meta"]:
                        mondo2desc[m_id] = n_obj["meta"]["definition"]["val"]


    newrow = ["xref_id","xref_key","name_type", "name", "description"]
    print "\"%s\""  % ("\",\"".join(newrow))
    for m_id in mondo2name:
        mondo2name[m_id] = mondo2name[m_id].encode('ascii', 'ignore').decode('ascii').replace("\n", "")
        desc = mondo2desc[m_id] if m_id in mondo2desc else ""
        desc = desc.encode('ascii', 'ignore').decode('ascii').replace("\"", "'").replace("\n", "")

        newrow = [m_id, "mondo", "recommended_name", mondo2name[m_id],desc]
        row_str = json.dumps(newrow)
        if row_str not in seen_row:
            print "\"%s\""  % ("\",\"".join(newrow))
            seen_row[row_str] = True
        if m_id in mondo2synonym:
            for syn in mondo2synonym[m_id]:
                syn = syn.encode('ascii', 'ignore').decode('ascii').replace("\n", "")
                newrow = [m_id, "mondo", "synonym", syn, ""]
                row_str = json.dumps(newrow)
                if row_str not in seen_row:
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen_row[row_str] = True


    data_grid = {
        "mimid2doid":{},
        "doid2icd10cm":{},
        "doid2icd9cm":{},
        "doid2kegg":{},
        "doid2mesh":{},
        "doid2umls":{},
        "doid2name":{},
        "doid2altname":{},
        "doid2def":{},
        "doid2mondo":{},
        "mimid2diseasename":{}
    }
    sparqlutil.load_do_mapping(data_grid)
    for species in species_list:
        sparqlutil.load_mimid2disease_name(data_grid, species)

    

    map_name = "mimid2diseasename"
    for k in data_grid[map_name]:
        desc = data_grid[map_name][k]
        desc = desc.encode('ascii', 'ignore').decode('ascii').replace("\"", "'")
        newrow = [k, "omim", "recommended_name", data_grid[map_name][k],desc]
        row_str = json.dumps(newrow)
        if row_str not in seen_row:
            print "\"%s\""  % ("\",\"".join(newrow).replace("\n", ""))
            seen_row[row_str] = True

    map_name = "doid2name"
    for do_id in data_grid[map_name]:
        for do_name in data_grid[map_name][do_id]:
            if do_name.strip() == "":
                continue
            do_name = do_name.encode('ascii', 'ignore').decode('ascii')
            desc = data_grid["doid2def"][do_id][0] if do_id in data_grid["doid2def"] else ""
            desc = desc.encode('ascii', 'ignore').decode('ascii').replace("\"", "'")
            newrow = [do_id, "do", "recommended_name", do_name,desc]
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                print "\"%s\""  % ("\",\"".join(newrow).replace("\n", ""))
                seen_row[row_str] = True

    map_name = "doid2altname"
    for do_id in data_grid[map_name]:
        for do_name in data_grid[map_name][do_id]:
            if do_name.strip() == "":
                continue
            do_name = do_name.encode('ascii', 'ignore').decode('ascii')
            newrow = [do_id, "do", "synonym", do_name, ""]
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                print "\"%s\""  % ("\",\"".join(newrow).replace("\n", ""))
                seen_row[row_str] = True


    return




def extract_disease_idmap_ds():

    #mondo_dict = {}
    seen_row = {}
    xrefid2mondoid = {}
    xrefid2doid = {}
    mondoid2doid = {}
    doid2mondoid = {}

    mondo_obj = json.loads(open(path_obj["downloads"] + "ohsu/mondo.json", "r").read())
    for g_obj in mondo_obj["graphs"]:
        for n_obj in g_obj["nodes"]:
            mondo_id = n_obj["id"].split("/")[-1]
            if mondo_id[0:6] == "MONDO_":
                m_id = mondo_id[6:]
                if "meta" in n_obj:
                    if "xrefs" in n_obj["meta"]:
                        for xref in n_obj["meta"]["xrefs"]:
                            xref_key = xref["val"].split(":")[0].lower()
                            xref_id = xref["val"].split(":")[-1]
                            #if xref_key not in ["OMIM", "DOID"]:
                            #    continue
                            #if m_id not in mondo_dict:
                            #    mondo_dict[m_id] = {}
                            #if xref_key not in mondo_dict[m_id]:
                            #    mondo_dict[m_id][xref_key] = []
                            #if xref_id not in  mondo_dict[m_id][xref_key]:
                            #    mondo_dict[m_id][xref_key].append(xref_id)
                            if xref_key not in xrefid2mondoid:
                                xrefid2mondoid[xref_key] = {}
                            if xref_id not in xrefid2mondoid[xref_key]:
                                xrefid2mondoid[xref_key][xref_id] = []
                            if m_id not in xrefid2mondoid[xref_key][xref_id]:
                                xrefid2mondoid[xref_key][xref_id].append(m_id)
                            if xref_key == "doid":
                                d_id = xref_id 
                                if m_id not in mondoid2doid:
                                    mondoid2doid[m_id] = []
                                if d_id not in mondoid2doid[m_id]:
                                    mondoid2doid[m_id].append(d_id)
                                
                                if d_id not in doid2mondoid:
                                    doid2mondoid[d_id] = []
                                if m_id not in doid2mondoid[d_id]:
                                    doid2mondoid[d_id].append(m_id)

    data_grid = { "mimid2doid":{}, "doid2icd10cm":{}, "doid2icd9cm":{}, "doid2kegg":{}, 
            "doid2mesh":{}, "doid2umls":{},
            "doid2name":{},"doid2altname":{},"doid2def":{}
    }
    
    sparqlutil.load_do_mapping(data_grid)
    xref_key = "omim"
    for xref_id in data_grid["mimid2doid"]:
        for do_id in data_grid["mimid2doid"][xref_id]:
            xref_id = xref_id.replace("PS", "")
            if xref_key not in xrefid2doid:
                xrefid2doid[xref_key] = {}
            if xref_id not in xrefid2doid[xref_key]:
                xrefid2doid[xref_key][xref_id] = []
            if do_id not in xrefid2doid[xref_key][xref_id]:
                xrefid2doid[xref_key][xref_id].append(do_id)

    row_list_one = []
    row_list_two = []
    for xref_id in xrefid2mondoid[xref_key]:
            m_list = xrefid2mondoid[xref_key][xref_id]
            d_list = xrefid2doid[xref_key][xref_id] if xref_id in xrefid2doid[xref_key] else [""]
            for m_id in m_list:
                if m_id in mondoid2doid and d_list == [""]:
                    d_list = mondoid2doid[m_id]
                    for d_id in d_list:
                        row_list_one.append([xref_id,xref_key,m_id,d_id,"mondo","mondo"])
                else:
                    for d_id in d_list:
                        if d_id != "":
                            row_list_one.append([xref_id,xref_key,m_id,d_id, "mondo", "do"])
                        else:
                            row_list_two.append([xref_id,xref_key,m_id,d_id, "mondo",""])


    row_list_three = []
    row_list_four = []
    for xref_id in xrefid2doid[xref_key]:
        d_list = xrefid2doid[xref_key][xref_id]
        m_list = xrefid2mondoid[xref_key][xref_id] if xref_id in xrefid2mondoid[xref_key] else [""]
        for d_id in d_list:
            if d_id in doid2mondoid and m_list == [""]:
                m_list = doid2mondoid[d_id]
                for m_id in m_list:
                    row_list_one.append([xref_id,xref_key,m_id,d_id,"mondo","do"])
            else:
                for m_id in m_list:
                    if m_id != "":
                        row_list_three.append([xref_id,xref_key,m_id,d_id,"mondo","do"])
                    else:
                        row_list_four.append([xref_id,xref_key,m_id,d_id,"mondo","do"])


    newrow = ["xref_id","xref_key","mondo_id", "do_id","xref_map_src", "do_map_src"]
    print "\"%s\""  % ("\",\"".join(newrow))

    is_complete = {}
    for row in row_list_one + row_list_three:
        row[1] = "protein_xref_" + row[1]
        combo_id = "%s|%s" % (row[1], row[0])
        is_complete[combo_id] = True
        row_str = json.dumps(row)
        if row_str not in seen_row:
            print "\"%s\""  % ("\",\"".join(row))
            seen_row[row_str] = True

    for row in row_list_two + row_list_four:
        row[1] = "protein_xref_" + row[1]
        combo_id = "%s|%s" % (row[1], row[0])
        if combo_id not in is_complete:
            row_str = json.dumps(row)
            if row_str not in seen_row:
                print "\"%s\""  % ("\",\"".join(row))
                seen_row[row_str] = True


    seen  = {}
    for in_file in glob.glob("compiled/*_protein_disease.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            xref_key  = row[f_list.index("xref_key")]
            xref_id  = row[f_list.index("xref_id")]
            if xref_key == "protein_xref_glygen_ds":
                do_id =  row[f_list.index("do_id")]
                if do_id not in seen:
                    if do_id not in doid2mondoid:
                        m_id = ""
                        newrow = [xref_id, xref_key, m_id, do_id,"mondo","do"]
                        print "\"%s\""  % ("\",\"".join(newrow))
                    else:
                        for m_id in doid2mondoid[do_id]:
                            newrow = [xref_id, xref_key, m_id, do_id, "mondo","do"]
                            print "\"%s\""  % ("\",\"".join(newrow))
                    seen[do_id] = True



    return



def extract_domap_dsold():


    data_grid = {
        "mimid2doid":{},
        "doid2icd10cm":{},
        "doid2icd9cm":{},
        "doid2kegg":{},
        "doid2mesh":{},
        "doid2umls":{},
        "doid2name":{},
        "doid2altname":{},
        "doid2def":{},
        "doid2mondo":{}
    }

    mondo_obj = json.loads(open(path_obj["downloads"] + "ohsu/mondo.json", "r").read())
    for g_obj in mondo_obj["graphs"]:
        for n_obj in g_obj["nodes"]:
            mondo_id = n_obj["id"].split("/")[-1]
            if "meta" in n_obj:
                if "xrefs" in n_obj["meta"]:
                    for xref in n_obj["meta"]["xrefs"]:
                        if mondo_id[0:6] == "MONDO_" and xref["val"][0:5] == "DOID:":
                            do_id = xref["val"][5:]
                            if do_id not in data_grid["doid2mondo"]:
                                data_grid["doid2mondo"][do_id] = []
                            data_grid["doid2mondo"][do_id].append(mondo_id[6:])
    sparqlutil.load_do_mapping(data_grid)

    row = ["do_id", "do_name", "xref_key", "xref_id"]
    print "\"%s\""  % ("\",\"".join(row))
            
    for do_id in data_grid["doid2name"]:
        for do_name in sorted(set(data_grid["doid2name"][do_id])):
            row = [do_id, do_name]
            for k in ["doid2kegg", "doid2icd10cm", "doid2icd9cm", "doid2mesh", "doid2umls", "doid2mondo"]:
                if do_id in data_grid[k]:
                    for db_id in sorted(set(data_grid[k][do_id])):
                        target = "protein_xref_" + k[5:]
                        newrow = row + [target, db_id]
                        print "\"%s\""  % ("\",\"".join(newrow))


    seen  = {}
    for in_file in glob.glob("compiled/*_protein_disease.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            xref_key  = row[f_list.index("xref_key")]
            xref_id  = row[f_list.index("xref_id")]
            if xref_key == "protein_xref_glygen_ds":
                do_id =  row[f_list.index("do_id")]
                if do_id not in seen:
                    do_name = ""
                    if do_id in data_grid["doid2name"]:
                        do_name = data_grid["doid2name"][do_id][0]
                    newrow = [do_id, do_name,xref_key, xref_id]
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen[do_id] = True

    seen = {}
    for mim_id in data_grid["mimid2doid"]:
        for do_id in data_grid["mimid2doid"][mim_id]:
            for do_name in sorted(set(data_grid["doid2name"][do_id])):
                mim_id = mim_id.replace("PS", "")
                combo_id = "%s %s" % (do_id, mim_id)
                if combo_id not in seen:
                    newrow = [do_id, do_name, "protein_xref_omim", mim_id]
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen[combo_id] = True





    return



def load_canon2xref(in_file, map_dict_one, map_dict_two):

    sheet_obj = {}
    libgly.load_sheet(sheet_obj, in_file, ",")
    f_list = sheet_obj["fields"]
    for row in sheet_obj["data"]:
        map_dict_one[row[f_list.index("xref_id")]] = row[f_list.index("uniprotkb_canonical_ac")]
        map_dict_two[row[f_list.index("uniprotkb_canonical_ac")]] = {
            "id":row[f_list.index("xref_id")],
            "name":row[f_list.index("xref_label")]
        }
    return



def extract_names_ds():
    lbl_dict = {
        "SemanticName":"Semantic Name",
        "ShortUniCarbKB":"UniCarbKB (Short)",
        "ShortComp":"Composition (Short)"
    }
    glycan_list = load_glycan_masterlist()
    in_file = path_obj["downloads"] + "glytoucan/current/export/names.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    newrow = ["glytoucan_ac","glycan_name","glycan_name_domain"]
    print "\"%s\""  % ("\",\"".join(newrow))
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        name = row[f_list.index("Name")]
        domain = row[f_list.index("Domain")]
        domain = lbl_dict[domain] if domain in lbl_dict else domain
        row[f_list.index("Domain")] = domain
        if ac in glycan_list:
            newrow = row
            print "\"%s\""  % ("\",\"".join(newrow))

    return

def extract_homolog_clusters_ds():


    work_book = {}
    xref2canon = {"geneid":{}, "hgnc":{}, "mgi":{}, "oma":{}}
    canon2xref = {"geneid":{}, "hgnc":{}, "mgi":{}, "oma":{}}

    
    xref = "hgnc"
    in_file = path_obj["unreviewed"] +  "human_protein_xref_hgnc.csv"
    load_canon2xref(in_file, xref2canon[xref], canon2xref[xref])

    xref = "mgi"
    in_file = path_obj["unreviewed"] +  "mouse_protein_xref_mgi.csv"
    load_canon2xref(in_file, xref2canon[xref], canon2xref[xref])

    xref = "oma"
    for in_file in glob.glob(path_obj["unreviewed"] + "*_protein_xref_oma.csv"):
        load_canon2xref(in_file, xref2canon[xref], canon2xref[xref])

    xref = "geneid"
    for in_file in glob.glob(path_obj["unreviewed"] + "*_protein_xref_geneid.csv"):
        load_canon2xref(in_file, xref2canon[xref], canon2xref[xref])

    
    pepid2canon = {}
    sheet_obj = {}
    for in_file in glob.glob(path_obj["unreviewed"] + "*_protein_transcriptlocus.csv"):
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        f_list = sheet_obj["fields"]
        for row in sheet_obj["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            uniprotkb_ac = canon.split("-")[0]
            isoform = row[f_list.index("uniprotkb_isoform_ac")]
            if isoform == canon:
                pep_id = row[f_list.index("peptide_id")]
                pepid2canon[pep_id] = canon
                pepid2canon[uniprotkb_ac] = canon


    out_rows = []
    seen_out_row = {}
    in_file = path_obj["downloads"] + "mgi/current/mgi_homologs.csv"
    if os.path.isfile(in_file) == True:
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        f_list = sheet_obj["fields"]
        homolog_dict = {}
        homologene2hgnc = {}
        homologene2mgi = {}
        homologene2taxid = {}
        for row in sheet_obj["data"]:
            #homologene_id = row[f_list.index("HomoloGene ID")]
            homologene_id = row[f_list.index("DB Class Key")]
            tax_id = row[f_list.index("NCBI Taxon ID")]
            gene_id = row[f_list.index("EntrezGene ID")]
            #hgnc_id = row[f_list.index("HGNC ID")].split(":")[-1]
            #mgi_id = row[f_list.index("Mouse MGI ID")].split(":")[-1]
            if gene_id in xref2canon["geneid"] and tax_id in ["9606", "10090", "10116"]:
                if homologene_id not in homologene2taxid:
                    homologene2taxid[homologene_id] = []
                if tax_id not in homologene2taxid[homologene_id]:
                    homologene2taxid[homologene_id].append(tax_id)
                canon = xref2canon["geneid"][gene_id]
                
                if homologene_id not in homolog_dict:
                    homolog_dict[homologene_id] = []
                homolog_dict[homologene_id].append({"canon":canon, "taxid": tax_id})

                #if tax_id == "9606":
                #    homologene2hgnc[homologene_id] = hgnc_id
                #elif tax_id == "10090":
                #    homologene2mgi[homologene_id] = mgi_id
        
        for homologene_id in homolog_dict:
            if homologene_id == "":
                continue
            #hgnc_id = homologene2hgnc[homologene_id] if homologene_id in homologene2hgnc else ""
            #mgi_id = homologene2mgi[homologene_id] if homologene_id in homologene2mgi else ""
            
            seen_taxid = {}
            homolog_set = []
            for obj in homolog_dict[homologene_id]:
                #if obj["taxid"] not in seen_taxid:
                #seen_taxid[obj["taxid"]] = True
                homolog_set.append(obj)
            
            if len(homolog_set) > 1:
                for obj in homolog_set:
                    tax_id = obj["taxid"]
                    tax_name = species_obj[tax_id]["long_name"]
                    row = [homologene_id, obj["canon"],tax_id, "protein_xref_mgi_homologset", 
                            homologene_id]
                    out_rows.append(row)

    for in_file in glob.glob(path_obj["downloads"] + "oma/current/*.csv"):
        species_one = in_file.split("/")[-1].split("_")[0]
        species_two = in_file.split("/")[-1].split("_")[1]
        tax_id_one = str(species_obj[species_one]["tax_id"])
        tax_id_two = str(species_obj[species_two]["tax_id"])
        tax_name_one = species_obj[species_one]["long_name"]
        tax_name_two = species_obj[species_two]["long_name"]

        peptide_map = {}
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        f_list = sheet_obj["fields"]
        for row in sheet_obj["data"]:

            oma_group = row[f_list.index("oma_group")]
            if oma_group == "":
                continue
            #if row[f_list.index("ortholog_subtype")] in ["1:1", "1:01"]:
            peptide_one = row[0].split(".")[0]
            peptide_two = row[1].split(".")[0]
            #print "Robel", row
            #print "Robel", peptide_one, peptide_two, peptide_one in pepid2canon, peptide_two in pepid2canon
            if peptide_one in pepid2canon and peptide_two in pepid2canon:
                canon_one, canon_two = pepid2canon[peptide_one], pepid2canon[peptide_two]
                out_rows.append([oma_group,canon_one,tax_id_one, "protein_xref_oma", oma_group])
                out_rows.append([oma_group,canon_two,tax_id_two,  "protein_xref_oma",oma_group])


    #GlyGen compiled clusters
    in_file = "compiled/protein_homolog_clusters.csv"
    libgly.load_sheet(sheet_obj, in_file, ",")
    f_list = sheet_obj["fields"]
    for row in sheet_obj["data"]:
        newrow = []
        for f in f_list:
            newrow.append(row[f_list.index(f)])
        out_rows.append(newrow)


    seen = {}
    cls_size = {}
    for row in out_rows:
        row_str = json.dumps(row)
        if row_str not in seen:
            seen[row_str]  = True
            cls = row[0]
            if cls not in cls_size:
                cls_size[cls] = 0
            cls_size[cls] += 1


    row = ["homolog_cluster_id", "uniprotkb_canonical_ac","tax_id", "xref_key", "xref_id"]
    print "\"%s\"" % ("\",\"".join(row))
    seen = {}
    for row in out_rows:
        cls = row[0]
        if cls_size[cls] < 2:
            continue
        row_str = json.dumps(row)
        if row_str not in seen:
            seen[row_str]  = True
            print "\"%s\"" % ("\",\"".join(row))



    return





def extract_homolog_alignments_ds():



    work_book = {}
    work_book["clusters"] = {}
    in_file = path_obj["unreviewed"] + "protein_homolog_clusters.csv"
    libgly.load_sheet(work_book["clusters"], in_file, ",")

    canon2cls = {}
    f_list = work_book["clusters"]["fields"]
    for row in work_book["clusters"]["data"]:
        database = row[f_list.index("xref_key")].split("_")[2]
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        cls = row[f_list.index("homolog_cluster_id")]
        cls = "%s.%s" % (database, cls)
        aln_file = "alignments/homologset/homologset.%s.%s.aln" % (database,cls)
        if os.path.isfile(aln_file) == True:
            file_size = int(os.path.getsize(aln_file))
            if file_size > 0:
                continue
        if canon not in canon2cls:
            canon2cls[canon] = []
        canon2cls[canon].append(cls)


    seq_hash = {}
    for fasta_file in glob.glob(path_obj["unreviewed"] +  "*protein_canon*sequences.fasta"):
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id.split("|")[1]
            seq_hash[seq_id] = str(record.seq.upper())

    seq_set = {}
    count_dict = {}
    for canon in canon2cls:
        if canon not in seq_hash:
            continue
        for cls in canon2cls[canon]:
            if cls not in seq_set:
                seq_set[cls] = ""
            seq_set[cls] += ">%s\n%s\n\n" % (canon, seq_hash[canon])
            if cls not in count_dict:
                count_dict[cls] = 0
            count_dict[cls] += 1


    ncls = 0
    for cls in seq_set:
        ncls += 1 if count_dict[cls] > 1 else 0

    with open("logs/homolog-cluster-aln.log", "w") as FL:
        FL.write("started making homolog cluster alignments\n")

    i = 0
    for cls in seq_set:
        if count_dict[cls] < 2:
            continue
        i += 1
        seq_file = "alignments/homologset/homologset.%s.fasta" % (cls)
        dnd_file = "alignments/homologset/homologset.%s.dnd" % (cls)
        aln_file = "alignments/homologset/homologset.%s.aln" % (cls)
        seq_buffer_old = ""
        if os.path.isfile(seq_file) == True:
            seq_buffer_old = open(seq_file, "r").read()
        seq_buffer_new = seq_set[cls] + "\n"
        if seq_buffer_new != seq_buffer_old or os.path.isfile(aln_file) == False:
            with open(seq_file, "w") as FW:
                FW.write("%s\n" % (seq_set[cls]))
            cmd = "/software/clustalw2/clustalw2 -align -type=protein -infile=%s " % (seq_file)
            x = commands.getoutput(cmd)
            cmd = "rm -f %s" % (dnd_file)
            x = commands.getoutput(cmd)
            with open("logs/homolog-cluster-aln.log", "a") as FL:
                FL.write("created %s\n" % (aln_file) )
        else:
            with open("logs/homolog-cluster-aln.log", "a") as FL:
                FL.write("file %s already exists\n" % (aln_file))

    return


def extract_species_customized_neuac_neugc_ds():

    seen = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "glycan_monosaccharide_composition.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        if row[f_list.index("NeuAc")] != "0" or row[f_list.index("NeuGc")] != "0":
            seen[ac] = True
    
    data_frame = {}
    in_file = path_obj["unreviewed"] + "glycan_species.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    newrow = f_list
    print "\"%s\"" % ("\",\"".join(newrow)) 

    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        if ac in seen:
            newrow = row
            print "\"%s\"" % ("\",\"".join(newrow))


 

    return



def extract_rdfdata_ds():
    src_file = path_obj["downloads"] + "glytoucan/current/export/glycandata.rdf.gz"
    dst_file = path_obj["unreviewed"] + "glycan_rdfdata.rdf.gz"
    
    cmd = "cp %s %s" % (src_file, dst_file)
    x = commands.getoutput(cmd)

    cmd = "readlink -f " + src_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")


    return

def extract_subsumption_ds():

    glycan_list = []
    glytoucan_type_dict = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        glytoucan_type = row[f_list.index("glytoucan_type")]
        if ac not in glycan_list:
            glycan_list.append(ac)
        glytoucan_type_dict[ac] = glytoucan_type


    in_file = path_obj["downloads"] + "glytoucan/current/export/subsumption.tsv"
    rename_dict = {
        "Relationship":"relationship",
        "RelatedAccession":"related_accession"
    }

    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    newrow = []
    for f in f_list:
        f_new = rename_dict[f] if f in rename_dict else f
        newrow.append(f_new)
    newrow[0] = "glytoucan_ac"
    newrow += ["glytoucan_type"]
    print "\"%s\"" % ("\",\"".join(newrow))


    for row in data_frame["data"]:
        newrow = row
        if newrow[0] in glycan_list:
            relationship = row[f_list.index("Relationship")]
            related_accession = row[f_list.index("RelatedAccession")]
            glytoucan_type = glytoucan_type_dict[related_accession] if related_accession in glytoucan_type_dict else ""
            print "\"%s\"" % ("\",\"".join(newrow + [glytoucan_type]))
    

    return



def extract_image_details_ds():
    glycan_list = load_glycan_masterlist()
    
    
    in_file = path_obj["downloads"] + "glytoucan/current/export/images.tsv"
    rename_dict = {
        "Image-Size":"image_size", 
        "Image-Notation":"image_notation",
        "Image-Style":"image_style",
        "Image-Format":"image_format"
    }

    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    newrow = []
    for f in f_list:
        f_new = rename_dict[f] if f in rename_dict else f
        newrow.append(f_new)
    newrow[0] = "glytoucan_ac"
    print "\"%s\"" % ("\",\"".join(newrow))

    for row in data_frame["data"]:
        newrow = row
        if newrow[0] in glycan_list:
            print "\"%s\"" % ("\",\"".join(newrow))

    return



def extract_images_ds(ds_name):



    reldir_dict = {
        "images_snfg_extended_png": "glytoucan/current/export/snfg/extended/png/"
        ,"images_snfg_extended_svg": "glytoucan/current/export/snfg/extended/svg/"
        ,"images_cfg_extended_png": "glytoucan/current/export/cfg/extended/png/"
        ,"images_cfg_extended_svg": "glytoucan/current/export/cfg/extended/svg/"
    }
    images_dir = path_obj["downloads"] + reldir_dict[ds_name]
    

    #Remove intermediate/glycanimages dir and create it again
    cmd = "rm -rf %s/%s" % (path_obj["intermediate"], ds_name)
    x = commands.getoutput(cmd)
    
    cmd = "mkdir %s/%s" % (path_obj["intermediate"], ds_name)
    x = commands.getoutput(cmd)

    #placeholder image
    cmd = "cp generated/misc/G00000*.png %s/%s/" % (path_obj["intermediate"], ds_name)
    x = commands.getoutput(cmd)
                

    seen_list = []
    data_frame = {}
    in_file = path_obj["unreviewed"] + "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        if ac not in seen_list:
            seen_list.append(ac)
            cmd = "cp %s/%s.* %s/%s/" % (images_dir, ac, path_obj["intermediate"], ds_name)
            x = commands.getoutput(cmd)           
    

    #Clean old dataset tarball
    cmd = "rm -f %s/glycan_%s*" % (path_obj["unreviewed"],ds_name)
    x = commands.getoutput(cmd)

    #Make tarball
    cmd = "/usr/bin/tar -C %s/%s/ -cvf %s/glycan_%s.tar ./" % (path_obj["intermediate"], ds_name, path_obj["unreviewed"], ds_name)
    x = commands.getoutput(cmd)
    cmd = "/usr/bin/gzip %s/glycan_%s.tar" % (path_obj["unreviewed"],ds_name)
    x = commands.getoutput(cmd)

    return

def extract_glytoucan_linkout_ds():

    newrow = ["id", "sequence"]
    print "\"%s\"" % ("\",\"".join(newrow))

    in_file = "unreviewed/glycan_sequences_wurcs.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for newrow in data_frame["data"]:
        print "\"%s\"" % ("\",\"".join(newrow))

    return


def get_ncbi_pubmed_linkouts_rows():

    row_list = []
    seen_row = {}
    for in_file in glob.glob("unreviewed/*_*_citations_*.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            if xref_key != "protein_xref_pubmed":
                continue
            url = "https://glygen.org/publication/PubMed/%s/" % (xref_id)
            newrow = ["10227", "PubMed",xref_id, url,"","","",""]
            if json.dumps(newrow) not in seen_row:
                row_list.append(newrow)
                seen_row[json.dumps(newrow)] = True

    return row_list


def extract_ncbi_pubmed_linkouts_one_ds():

    row_list = get_ncbi_pubmed_linkouts_rows()
    batch_size = int(len(row_list)/3)
    start = 0*batch_size
    end = start + batch_size
    newrow = ["ProviderId","Database","UID", "URL", "IconUrl", "UrlName", 
            "SubjectType", "Attribute"]
    print "%s" % (",".join(newrow))
    for newrow in row_list[start:end]:
        print "%s" % (",".join(newrow))

    return

def extract_ncbi_pubmed_linkouts_two_ds():

    row_list = get_ncbi_pubmed_linkouts_rows()
    batch_size = int(len(row_list)/3)
    start = 1*batch_size
    end = start + batch_size
    newrow = ["ProviderId","Database","UID", "URL", "IconUrl", "UrlName", 
            "SubjectType", "Attribute"]
    print "%s" % (",".join(newrow))
    for newrow in row_list[start:end]:
        print "%s" % (",".join(newrow))
    
    return


def extract_ncbi_pubmed_linkouts_three_ds():

    row_list = get_ncbi_pubmed_linkouts_rows()
    batch_size = int(len(row_list)/3)
    start = 2*batch_size
    newrow = ["ProviderId","Database","UID", "URL", "IconUrl", "UrlName", 
            "SubjectType", "Attribute"]
    print "%s" % (",".join(newrow))
    for newrow in row_list[start:]:
        print "%s" % (",".join(newrow))
    
    return



def extract_ncbi_gene_linkouts_ds():

    newrow = ["ProviderId","Database","UID", "URL", "IconUrl", "UrlName", "SubjectType", "Attribute"]
    print "%s" % (",".join(newrow))
    seen_row = {}
    count_dict = {}
    for in_file in glob.glob("unreviewed/*_protein_xref_geneid.csv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            uniprotkb_ac = row[f_list.index("uniprotkb_canonical_ac")].split("-")[0]
            gene_id = row[f_list.index("xref_id")]
            if gene_id not in count_dict:
                count_dict[gene_id] = 0
            count_dict[gene_id] += 1
            if count_dict[gene_id] < 4:
                url = "https://glygen.org/protein/%s" % (uniprotkb_ac)
                newrow = ["10227", "Gene",gene_id, url,"","","",""]
                if json.dumps(newrow) not in seen_row:
                    print "%s" % (",".join(newrow))
                    seen_row[json.dumps(newrow)] = True 

    return 


def extract_uniprotkb_accession_history_ds():

    newrow = ["uniprotkb_ac_old", "uniprotkb_ac_current"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for in_file in glob.glob(path_obj["downloads"] + "ebi/current/accession-history-*.tsv"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, "\t")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            old_ac = row[f_list.index("old_accession")]
            current_ac = row[f_list.index("current_accession")] if len(row) > 1  else ""
            for ac_one in old_ac.strip().split(","):
                for ac_two in current_ac.strip().split(","):
                    if ac_one.strip() != "" and ac_two.strip() != "":
                        newrow = [ac_one.strip(), ac_two.strip()]
                        print "\"%s\"" % ("\",\"".join(newrow))


    return


def extract_pubchem_status_ds():

    newrow = ["glytoucan_ac","pubchem_xref_exists", "xref_key", "xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    data_frame = {}
    in_file = path_obj["unreviewed"] + "glycan_xref_pubchem.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    status_dict = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")]   
        if glytoucan_ac not in status_dict:
            status_dict[glytoucan_ac] = {}
        combo = row[f_list.index("xref_key")] + "|"  + row[f_list.index("xref_id")] 
        status_dict[glytoucan_ac][combo] =  True
    
    data_frame = {}
    in_file = path_obj["unreviewed"] + "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")]
        if glytoucan_ac in status_dict:
            xref_key_list, xref_id_list = [], []
            for combo in status_dict[glytoucan_ac]:
                xref_key, xref_id = combo.split("|")
                xref_key_list.append(xref_key)
                xref_id_list.append(xref_id)
            newrow = [glytoucan_ac, "yes", "|".join(xref_key_list), "|".join(xref_id_list)]
            print "\"%s\"" % ("\",\"".join(newrow))
        else:
            newrow = [glytoucan_ac, "no", "", ""]
            print "\"%s\"" % ("\",\"".join(newrow))
            

    return



def extract_pathway_reactome_ds():


    newrow = ["glytoucan_ac","compound_name", "reactome_pathway_id","pathway_name",
            "species_scientific_name", "xref_key", "xref_id"] 
    print "\"%s\"" % ("\",\"".join(newrow))

    cid2gtc = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "glycan_xref_chebi.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")]
        c_id = row[f_list.index("xref_id")]
        if c_id not in cid2gtc:
            cid2gtc[c_id] = []
        if glytoucan_ac not in cid2gtc[c_id]:
            cid2gtc[c_id].append(glytoucan_ac)

    data_frame = {}
    in_file = path_obj["downloads"] + "reactome/current/ChEBI2Reactome_PE_Pathway.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    k_list = ["chebi_id","reactome_compound_id", "compound_name", 
            "reactome_pathway_id","pathway_name", "species_scientific_name"]
    for row in data_frame["data"]:
        val_dict = {}
        for k in k_list:
            val_dict[k] = row[f_list.index(k)]
        c_id = val_dict["chebi_id"]
        if c_id not in cid2gtc:
            continue
        xref_key, xref_id = "glycan_xref_reactome", val_dict["reactome_compound_id"]
        for glytoucan_ac in cid2gtc[c_id]:
            newrow = [glytoucan_ac]
            for k in k_list[2:]:
                newrow.append(val_dict[k])
            newrow += [xref_key, xref_id]
            print "\"%s\"" % ("\",\"".join(newrow))
    return








###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-d","--dataset",action="store",dest="dataset",help="[masterlist, transcriptlocus]")

    (options,args) = parser.parse_args()
    for file in ([options.dataset]):
        if not (file):
            parser.print_help()
            sys.exit(0)


    start_time = time.time()


    global config_obj
    global sparql
    global graph_uri
    global prefixes
    global data_grid
    global species_obj


    dataset = options.dataset

    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = {}
    in_file = config_obj["pathinfo"]["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)



    global path_obj
    path_obj = config_obj["pathinfo"]


    
    libgly.log_file_usage("", dataset, "write")


    mol = "glycan"
    #glycan datasets
    if dataset == "glytoucanidlist":
        extract_glytoucanidlist_ds()
    elif dataset == "masterlist":
        extract_masterlist_ds()
    elif dataset == "xref_chebi":
        extract_xref_chebi_ds()
    elif dataset == "enzyme":
        extract_enzyme_ds()
    elif dataset == "species":
        extract_species_ds()
    elif dataset == "fully_determined":
        extract_fully_determined_ds()
    elif dataset == "citations_glytoucan":
        extract_citations_glytoucan_ds()
    elif dataset == "citations_ncfg":
        extract_citations_ncfg_ds()
    elif dataset == "citations_glycomotif":
        extract_citations_glycomotif_ds()
    elif dataset == "classification":
        extract_classification_ds()
    elif dataset == "motif":
        extract_motif_ds()
    elif dataset == "monosaccharide_composition":
        extract_monosaccharide_composition_ds()
    elif dataset == "monosaccharide_composition_advanced":
        extract_monosaccharide_composition_advanced_ds()
    elif dataset in ["sequences_iupac_extended", "sequences_glycoct", "sequences_wurcs", 
            "sequences_glycam_iupac", "sequences_byonic", "sequences_gwb"]:
        extract_sequences_ds("_".join(dataset.split("_")[1:]) )
    elif dataset == "sequences_smiles_isomeric":
        extract_sequences_smiles_isomeric_ds()
    elif dataset == "sequences_inchi":
        extract_sequences_inchi_ds()
    elif dataset == "names":
        extract_names_ds()
    elif dataset in config_obj["xref"]:
        extract_xrefs_ds(dataset)
    elif dataset == "pubchem_status":
        extract_pubchem_status_ds()
    elif dataset in ["images_snfg_extended_png", "images_snfg_extended_svg"]:
        extract_images_ds(dataset)
    elif dataset == "image_details":
        extract_image_details_ds()
    elif dataset == "evidence_ncfg":
        extract_evidence_ncfg_ds()
    elif dataset in  ["customized_neuac_mammalian_species_annotation", "synthesized"]:
        extract_compiled_ds(dataset, True)
    elif dataset == "subsumption":
        extract_subsumption_ds()
    elif dataset in  ["rdfdata"]:
        extract_rdfdata_ds()
    elif dataset == "species_customized_neuac_neugc":
        extract_species_customized_neuac_neugc_ds()
    elif dataset == "type_n_linked_byonic":
        extract_type_n_linked_byonic_ds()
    elif dataset == "glytoucan_linkout":
        extract_glytoucan_linkout_ds()
    elif dataset == "glytoucan_accession_history":
        extract_glytoucan_accession_history_ds()
    elif dataset == "pathway_reactome":
        extract_pathway_reactome_ds()
    elif dataset == "pathway_glycotree":
        extract_pathway_glycotree_ds()
    elif dataset == "top_authors":
        extract_top_authors_ds()
    elif dataset == "dictionary":
        extract_dictionary_ds()




    #protein datasets
    if dataset == "homolog_clusters":
        mol = "protein"
        extract_homolog_clusters_ds()
    elif dataset == "homolog_alignments":
        mol = "protein"
        extract_homolog_alignments_ds()
    elif dataset == "disease_idmap":
        mol = "protein"
        extract_disease_idmap_ds()
    elif dataset == "disease_names":
        mol = "protein"
        extract_disease_names_ds()
    elif dataset == "ncbi_gene_linkouts":
        mol = "protein"
        extract_ncbi_gene_linkouts_ds()
    elif dataset == "ncbi_pubmed_linkouts_one":
        mol = "protein"
        extract_ncbi_pubmed_linkouts_one_ds()
    elif dataset == "ncbi_pubmed_linkouts_two":
        mol = "protein"
        extract_ncbi_pubmed_linkouts_two_ds()
    elif dataset == "ncbi_pubmed_linkouts_three":
        mol = "protein"
        extract_ncbi_pubmed_linkouts_three_ds()
    elif dataset in ["glygen_pharos_xref_mapping", "glygen_pubchem_xref_mapping",
            "glygen_uniprotkb_xref_mapping", "glygen_genecards_xref_mapping", 
            "glygen_iptmnet_xref_mapping"]:
        mol = "protein"
        extract_xref_mapping_ds(dataset)
    elif dataset == "glygen_uniprotkb_protvista_mapping":
        mol = "protein"
        extract_glygen_uniprotkb_protvista_mapping_ds()
    elif dataset == "uniprotkb_accession_history":
        mol = "protein"
        extract_uniprotkb_accession_history_ds()


    pid = os.getpid()
    pid = os.getpid()
    src_file = "usage/file_usage.%s.log" % (pid)
    dst_file = "usage/%s_%s.fu.log" % (mol, dataset)

    cmd = "cat %s |sort -u > %s" % (src_file, dst_file)
    x = commands.getoutput(cmd)
    
    cmd = "rm -f %s" % (src_file)
    x = commands.getoutput(cmd)
    



if __name__ == '__main__':
        main()

