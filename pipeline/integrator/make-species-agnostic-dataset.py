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
#import requests
import gzip


sys.path.append('../../glytools/')
import libgly


__version__="1.0"
__status__ = "Dev"




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

    newrow = ["glytoucan_ac","sequence_%s" % (seq_format)]
    print "\"%s\"" % ("\",\"".join(newrow))

    seen_list = []
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        if ac not in seen_list:
            seen_list.append(ac)
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






def extract_motif_ds():

    glycan_list = load_glycan_masterlist()
    newrow = ["glytoucan_ac","glytoucan_ac_motif","motif_name", "is_reducing_end"]
    print "\"%s\"" % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/allmotif.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    seen = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        motif_ac = row[f_list.index("MotifAccession")]
        motif_label = row[f_list.index("Label")]
        is_reducing_end = row[f_list.index("IsReducingEnd")]
        combo_id = "%s %s" % (ac, motif_ac)
        if ac in glycan_list:
            if combo_id not in seen:
                seen[combo_id] = True
                newrow = [ac, motif_ac, motif_label, is_reducing_end]
                print "\"%s\"" % ("\",\"".join(newrow))


    return



def extract_classification_ds():

    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","glycan_type","glycan_subtype"]
    print "\"%s\"" % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/classification.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        if glytoucan_ac in glycan_list:
            glycan_type = row[f_list.index("Type")]
            glycan_type = "N-glycan" if glycan_type == "N-linked" else glycan_type
            glycan_type = "O-glycan" if glycan_type == "O-linked" else glycan_type
            glycan_subtype = row[f_list.index("Subtype")]
            newrow = [glytoucan_ac, glycan_type, glycan_subtype]
            print "\"%s\"" % ("\",\"".join(newrow))


    return



def extract_citations_glytoucan_ds():

    black_list = get_blacklisted_pmids()

    glycan_list = load_glycan_masterlist()

    in_file = path_obj["downloads"] + "glytoucan/current/export/pubs.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    newrow = ["glytoucan_ac","pmid","title","journal_name","publication_date", "authors","source", "source_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        pmid = row[f_list.index("PubMedID")]
        source = row[f_list.index("Source")]
        source_id = row[f_list.index("SourceID")]
        if source_id.find("comp_") != -1:
            continue
        if pmid in black_list:
            continue
        combo_id = "%s %s" % (glytoucan_ac, pmid) 
        cond_list = [glytoucan_ac not in glycan_list and glytoucan_ac != ""]
        cond_list.append(pmid in ["0"])
        cond_list.append(combo_id in seen)
        if True in cond_list:
            continue
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
            if "TI" not in obj or "JT" not in obj:
                continue
            title = " ".join(obj["TI"])
            journal = " ".join(obj["JT"])
            pubdate = " ".join(obj["DP"])
            authors = ", ".join(obj["AU"])
            newrow = [glytoucan_ac, pmid, title, journal, pubdate, authors, source, source_id]
            print "\"%s\"" % ("\",\"".join(newrow))
            seen[combo_id] = True


    return




def extract_masterlist_ds():

    seen_list = []
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/allglycan.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        if ac not in seen_list:
            seen_list.append(ac)

    newrow = [
        "glytoucan_ac","glytoucan_type","glycan_mass", "glycan_permass",
        "base_composition","composition","topology","monosaccharides"
    ]
    print "\"%s\"" % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/glycan_properties.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
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
        if newrow[0] in seen_list:
            print "\"%s\"" % ("\",\"".join(newrow))

    return



def extract_fully_determined_ds():

    glycan_list = load_glycan_masterlist()
    newrow = ["glytoucan_ac"]
    print "\"%s\""  % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/fully_determined.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        if ac in glycan_list:
            newrow = [ac]
            print "\"%s\""  % ("\",\"".join(newrow))

    return




def extract_taxonomy_ds():


    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","tax_id", "source", "source_id"]
    print "\"%s\""  % ("\",\"".join(newrow))

    seen = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/taxa.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        tax_id = row[f_list.index("TaxID")]
        source = row[f_list.index("Source")]
        source_id = row[f_list.index("SourceID")]
        if source_id.find("comp_") != -1:
            continue
        if ac in glycan_list:
            newrow = [ac, tax_id, source, source_id]
            newrow_str = ",".join(newrow)
            if newrow_str not in seen:
                print "\"%s\""  % ("\",\"".join(newrow))
                seen[newrow_str] = True

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/species.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("GlyTouCanAccession")]
        evidence_type = row[f_list.index("Species")]
        evidence_desc = row[f_list.index("Value")].replace("[", "").replace("]", "").replace("'", "")
        if evidence_desc.find("TaxID") == -1:
            continue
        newrow = []
        if evidence_desc.find("GlyTouCan") != -1:
            tax_id = evidence_desc.split(" ")[-1]
            source = "GlyTouCan"
            source_id = ac
            newrow = [ac,tax_id, source, source_id]
        elif evidence_desc.find("UniCarbKB") != -1:
            tax_id = evidence_desc.split(" ")[-1]
            source = "UniCarbKB"
            source_id = evidence_desc.split(" ")[-3].split(":")[1]
            if source_id.find("comp_") == -1:
                newrow = [ac,tax_id, source, source_id]
        if ac in glycan_list and newrow != []:
            newrow_str = ",".join(newrow)
            if newrow_str not in seen:
                print "\"%s\""  % ("\",\"".join(newrow))
                seen[newrow_str] = True


    return


def extract_synthesized_ds():

    glycan_list = load_glycan_masterlist()

    in_file = "compiled/glycan_synthesized.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    newrow = f_list
    print "\"%s\""  % ("\",\"".join(newrow))
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        if ac in glycan_list:
            newrow = row
            print "\"%s\""  % ("\",\"".join(newrow))

    return


def extract_xrefs_ds(dataset):


    glycan_list = load_glycan_masterlist()

    newrow = ["glytoucan_ac","database_id","database_label"]
    print "\"%s\""  % ("\",\"".join(newrow))

    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/%s" % (config_obj["xref"][dataset]["glycaninfile"])
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[0]
        col_ind = config_obj["xref"][dataset]["colind"] if "colind" in config_obj["xref"][dataset] else 1
        database_id = row[col_ind]
        if database_id == "":
            continue
        if ac in glycan_list:
            db_name = "KEGG Glycan" if dataset == "xref_kegg" else config_obj["xref"][dataset]["dbname"]
            newrow = [ac, database_id, db_name]
            print "\"%s\""  % ("\",\"".join(newrow))

    return




def extract_enzyme_ds():


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


    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/enz.tsv"
    if os.path.isfile(in_file) == False:
        return
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    FL = open(path_obj["logs"] + "glycan_enzyme.csv" , "w")
    enz_dict = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        uniprotkb_ac = row[f_list.index("Enzyme_UniProtAccession")]
        gene_id = row[f_list.index("Enzyme_GeneID")]
        tax_id = row[f_list.index("Enzyme_TaxonomyID")]
        if glytoucan_ac not in enz_dict:
            enz_dict[glytoucan_ac] = []
        if uniprotkb_ac not in ac2canon:
            FL.write("\"flag-1\",\"%s\"\n" % ("\",\"".join(row)))
        else: 
            canon = ac2canon[uniprotkb_ac]
            gene_name = canon2genename[canon] if canon in canon2genename else ""
            recname = canon2recname[canon] if canon in canon2recname else ""
            o = {
                "canon":canon,
                "recname":recname,
                "gene_name":gene_name,
                "glytoucan_ac":glytoucan_ac,
                "gene_id":gene_id,
                "tax_id":tax_id
            }
            enz_dict[glytoucan_ac].append(o)
    

    #[GlyTouCanAccession,Enzyme_UniProtAccession, Enzyme_GeneName,Enzyme_GeneID,Enzyme_TaxonomyID]
    #[glycan_ID,residue,residue_ID,parent_ID]
    #[inputfile1: GlyTouCanAccession ; inputfile2: glycan_ID]


    data_frame = {}
    in_file = path_obj["downloads"] + "uga/mapped_n_glycans.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    newrow = ["glytoucan_ac","residue_name","residue_id","uniprotkb_canonical_ac","gene_name","gene_id","parent_residue_id","recommended_name_full"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glycan_ID")]
        res_name = row[f_list.index("residue")]
        res_id = row[f_list.index("residue_ID")]
        parent_id = row[f_list.index("parent_ID")]
        if glytoucan_ac in enz_dict:
            for o in enz_dict[glytoucan_ac]:
                newrow = [glytoucan_ac,res_name,res_id,o["canon"],o["gene_name"],o["gene_id"],parent_id,o["recname"]]
                print "\"%s\"" % ("\",\"".join(newrow))
    FL.close()



def extract_xref_chebi_from_kegg_ds():


    kegg2chebi = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "chebi/database_accession_current.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        if row[f_list.index("SOURCE")] == "KEGG GLYCAN":
            chebi_id = row[f_list.index("COMPOUND_ID")]
            kegg_id = row[f_list.index("ACCESSION_NUMBER")]
            kegg2chebi[kegg_id] = chebi_id


    newrow = ["glytoucan_ac","database_id","database_label"]
    print "\"%s\"" % ("\",\"".join(newrow))

    in_file = path_obj["unreviewed"] +  "glycan_xref_kegg.csv" 
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        kegg_id = row[f_list.index("database_id")]
        if kegg_id in kegg2chebi:
            chebi_id = kegg2chebi[kegg_id]
            newrow = [ac,chebi_id,"ChEBI"]
            print "\"%s\"" % ("\",\"".join(newrow))



    return





def extract_xref_chebi_ds():

    chebi2inchi = {}
    in_file = path_obj["downloads"] + "chebi/chebiId_inchi.tsv"
    with open(in_file, "r") as FR:
        for line in FR:
            parts = line.strip().split("\t")
            chebi2inchi[parts[0]] = parts[1]
    
    pubchem2inchi = {}
    in_file = "generated/pubchem/compound/cid2inchi.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        pubchem_id = row[f_list.index("pubchem_cid")]
        inchi = row[f_list.index("inchi")]
        pubchem2inchi[pubchem_id] = inchi
    

    cid2glytoucan = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_xref_pubchem.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        database_id = row[f_list.index("database_id")]
        if database_id[0:3] == "CID":
            cid2glytoucan[database_id[3:]] = ac

    newrow = ["glytoucan_ac","database_id","database_label"]
    print "\"%s\""  % ("\",\"".join(newrow))
    in_file = path_obj["downloads"] + "pubchem/compound/cid2synonym.tsv"
    with open(in_file, "r") as FR:
        for line in FR:
            line = line.strip()
            cid, chebi_id = line.split("\t")
            if cid in cid2glytoucan and chebi_id[0:6] == "CHEBI:":
                if chebi_id[6:] in chebi2inchi and cid in pubchem2inchi:
                    chebi_inchi = chebi2inchi[chebi_id[6:]]
                    pubchem_inchi = pubchem2inchi[cid]
                    #newrow = [cid2glytoucan[cid], "CID" + cid, chebi_id]
                    newrow = [cid2glytoucan[cid], chebi_id.split(":")[1], "ChEBI"]
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
        database_id = row[f_list.index("database_id")]
        if database_id[0:3] == "CID":
            cid2glytoucan[database_id[3:]] = ac

    newrow = ["glytoucan_ac","pubchem_id", "sequence_smiles_isomeric"]
    print "\"%s\"" % ("\",\"".join(newrow))
    in_file = path_obj["downloads"] + "pubchem/compound/cid2smiles.tsv"
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
        database_id = row[f_list.index("database_id")]
        if database_id[0:3] == "CID":
            cid2glytoucan[database_id[3:]] = ac

    newrow = ["glytoucan_ac","sequence_inchi","inchi_key"]
    print "\"%s\"" % ("\",\"".join(newrow))

    data_frame = {}
    in_file = "generated/pubchem/compound/cid2inchi.csv"
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



def extract_domap_ds():


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
    
    row = ["do_id", "do_name", "xref_database", "xref_database_id"]
    print "\"%s\""  % ("\",\"".join(row))
            
    for do_id in data_grid["doid2name"]:
        for do_name in sorted(set(data_grid["doid2name"][do_id])):
            row = [do_id, do_name]
            for k in ["doid2kegg", "doid2icd10cm", "doid2icd9cm", "doid2mesh", "doid2umls", "doid2mondo"]:
                if do_id in data_grid[k]:
                    for db_id in sorted(set(data_grid[k][do_id])):
                        target = k[5:]
                        newrow = row + [target, db_id]
                        print "\"%s\""  % ("\",\"".join(newrow))


    for mim_id in data_grid["mimid2doid"]:
        for do_id in data_grid["mimid2doid"][mim_id]:
            for do_name in sorted(set(data_grid["doid2name"][do_id])):
                mim_id = mim_id.replace("PS", "")
                newrow = [do_id, do_name, "omim", mim_id]
                print "\"%s\""  % ("\",\"".join(newrow))


    return



def load_canon2xref(in_file, map_dict_one, map_dict_two):

    sheet_obj = {}
    libgly.load_sheet(sheet_obj, in_file, ",")
    f_list = sheet_obj["fields"]
    for row in sheet_obj["data"]:
        map_dict_one[row[f_list.index("database_id")]] = row[f_list.index("uniprotkb_canonical_ac")]
        map_dict_two[row[f_list.index("uniprotkb_canonical_ac")]] = {
            "id":row[f_list.index("database_id")],
            "name":row[f_list.index("database_label")]
        }
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
            isoform = row[f_list.index("uniprotkb_isoform_ac")]
            if isoform == canon:
                pep_id = row[f_list.index("peptide_id")]
                pepid2canon[pep_id] = canon

    out_rows = []
    seen_out_row = {}
    in_file = path_obj["downloads"] + "mgi/current.csv"
    if os.path.isfile(in_file) == True:
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        f_list = sheet_obj["fields"]
        homolog_dict = {}
        homologene2hgnc = {}
        homologene2mgi = {}
        homologene2taxid = {}
        for row in sheet_obj["data"]:
            homologene_id = row[f_list.index("HomoloGene ID")]
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
                    tax_name = species_obj[tax_id]["taxname"]
                    row = [homologene_id, obj["canon"],tax_id, tax_name, "mgi"]
                    out_rows.append(row)

    for in_file in glob.glob(path_obj["downloads"] + "oma/current/*.csv"):
        species_one = in_file.split("/")[-1].split("_")[0]
        species_two = in_file.split("/")[-1].split("_")[1]
        tax_id_one = str(species_obj[species_one]["taxid"])
        tax_id_two = str(species_obj[species_two]["taxid"])
        tax_name_one = species_obj[species_one]["taxname"]
        tax_name_two = species_obj[species_two]["taxname"]

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
            if peptide_one in pepid2canon and peptide_two in pepid2canon:
                canon_one, canon_two = pepid2canon[peptide_one], pepid2canon[peptide_two]
                out_rows.append([oma_group,canon_one,tax_id_one, tax_name_one, "oma"])
                out_rows.append([oma_group,canon_two,tax_id_two, tax_name_two, "oma"])

    row = ["homolog_cluster_id", "uniprotkb_canonical_ac","tax_id", "tax_name", "database"]
    print "\"%s\"" % ("\",\"".join(row))
    seen = {}
    for row in out_rows:
        row_str = json.dumps(row)
        if row_str not in seen:
            seen[row_str]  = True
            print "\"%s\"" % ("\",\"".join(row))



    return





def extract_homolog_alignments_ds():



    work_book = {}
    work_book["clusters"] = {}
    in_file = path_obj["unreviewed"] + "homolog_clusters.csv"
    libgly.load_sheet(work_book["clusters"], in_file, ",")

    canon2cls = {}
    f_list = work_book["clusters"]["fields"]
    for row in work_book["clusters"]["data"]:
        database = row[f_list.index("database")]
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        cls = row[f_list.index("homolog_cluster_id")]
        if database in ["oma", "mgi"]:
            cls = "%s_group_%s" % (database, cls)
            aln_file = "alignments/homologset.%s.%s.aln" % (database,cls)
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
        for cls in canon2cls[canon]:
            if cls not in seq_set:
                seq_set[cls] = ""
            seq_set[cls] += ">%s\n%s\n\n" % (canon, seq_hash[canon])
            if cls not in count_dict:
                count_dict[cls] = 0
            count_dict[cls] += 1




    for cls in seq_set:
        if count_dict[cls] < 2:
            continue
        seq_file = "alignments/homologset.%s.fasta" % (cls)
        dnd_file = "alignments/homologset.%s.dnd" % (cls)
        with open(seq_file, "w") as FW:
            FW.write("%s\n" % (seq_set[cls]))
        cmd = "/software/clustalw2/clustalw2 -align -type=protein -infile=%s " % (seq_file)
        x = commands.getoutput(cmd)
        cmd = "rm -f %s %s" % (seq_file, dnd_file)
        x = commands.getoutput(cmd)
        print "aligned cls=%s" % (cls)

    return


def extract_images_ds():

    images_dir = path_obj["downloads"] +  "glytoucan/current/export/images-cfg-extended/cfg/extended/"
        
    #Remove intermediate/glycanimages dir and create it again
    cmd = "rm -rf %s/glycanimages" % (path_obj["intermediate"])
    x = commands.getoutput(cmd)
    
    cmd = "mkdir %s/glycanimages" % (path_obj["intermediate"])
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
            cmd = "cp %s/%s.png %s/glycanimages/" % (images_dir, ac, path_obj["intermediate"])
            x = commands.getoutput(cmd)           

    cmd = "/usr/bin/tar -C %s/glycanimages/ -cvf %s/glycan_images.tar ./" % (path_obj["intermediate"], path_obj["unreviewed"])
    x = commands.getoutput(cmd)
    cmd = "/usr/bin/gzip %s/glycan_images.tar" % (path_obj["unreviewed"])
    x = commands.getoutput(cmd)

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


    global config_obj
    global sparql
    global graph_uri
    global prefixes
    global data_grid
    global species_obj


    dataset = options.dataset

    config_obj = json.loads(open("conf/config.json", "r").read())
    species_obj = config_obj["speciesinfo"]

    global path_obj
    path_obj = config_obj["pathinfo"]


    data_grid = {}
    if dataset == "masterlist":
        extract_masterlist_ds()
    elif dataset == "xref_chebi":
        extract_xref_chebi_ds()
    elif dataset == "xref_chebi_from_kegg":
        extract_xref_chebi_from_kegg_ds()
    elif dataset == "enzyme":
        extract_enzyme_ds()
    elif dataset == "taxonomy":
        extract_taxonomy_ds()
    elif dataset == "fully_determined":
        extract_fully_determined_ds()
    elif dataset == "citations_glytoucan":
        extract_citations_glytoucan_ds()
    elif dataset == "classification":
        extract_classification_ds()
    elif dataset == "motif":
        extract_motif_ds()
    elif dataset == "monosaccharide_composition":
        extract_monosaccharide_composition_ds()
    elif dataset in ["sequences_iupac_extended", "sequences_glycoct", "sequences_wurcs", 
            "sequences_glycam_iupac"]:
        extract_sequences_ds("_".join(dataset.split("_")[1:]) )
    elif dataset == "sequences_smiles_isomeric":
        extract_sequences_smiles_isomeric_ds()
    elif dataset == "sequences_inchi":
        extract_sequences_inchi_ds()
    elif dataset in config_obj["xref"]:
        extract_xrefs_ds(dataset)
    elif dataset == "homolog_clusters":
        extract_homolog_clusters_ds()
    elif dataset == "homolog_alignments":
        extract_homolog_alignments_ds()
    elif dataset == "domap":
        extract_domap_ds()
    elif dataset == "images":
        extract_images_ds()
    elif dataset == "synthesized":
        extract_synthesized_ds()






if __name__ == '__main__':
        main()

