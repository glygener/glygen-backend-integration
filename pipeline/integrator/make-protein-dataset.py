import os,sys
import json
import csv
#import requests

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


def expand_rowcoll(rowcoll, val_list):

    newrowcoll = []
    for row in rowcoll:
        for val in val_list:
            val = val.encode('ascii', 'ignore').decode('ascii') if val != None else ""
            newrowcoll.append(row + [val])
    return newrowcoll


def extract_genomics_england_disease_ds(species):

    row_list = get_genomics_england_disease_rowlist(species)
    for newrow in row_list:
        print "\"%s\"" % ("\",\"".join(newrow))
    return


def get_genomics_england_disease_rowlist(species):

    mimid2doid = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "protein_domap.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        xref_database = row[f_list.index("xref_database")]
        xref_database_id = row[f_list.index("xref_database_id")]
        if xref_database == "omim":
            mimid2doid[xref_database_id] = do_id


    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_hgnc.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    genename2canon = {}
    genename2hgncid = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_name = row[f_list.index("database_label")]
        genename2canon[gene_name] = canon
        genename2hgncid[gene_name] = row[f_list.index("database_id")]


    row_list = [
        ["uniprotkb_canonical_ac","gene_name", "hgnc_id", "mode_of_inheritance","disease_name",
            "pmid","mim_id","doid"]
    ]

    data_frame = {}
    in_file = "downloads/genomics_england/current/congenital_disorders_of_glycosylation.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        gene_name = row[f_list.index("Gene Symbol")]
        inheritance = row[f_list.index("Model_Of_Inheritance")]
        phenotype_list = row[f_list.index("Phenotypes")].split(";")
        pmid_list = row[f_list.index("Publications")].split(";")
        sources = row[f_list.index("Sources(; separated)")]

        if sources.find("Expert Review Green") == -1:
            continue
        if gene_name not in genename2canon:
            continue
        if gene_name not in genename2hgncid:
            continue
        canon = genename2canon[gene_name]
        hgnc_id = genename2hgncid[gene_name]
        for p in phenotype_list:
            disease_name = " ".join(p.split(" ")[:-1])
            mim_id = p.split(" ")[-1]
            if mim_id in mimid2doid:
                do_id = mimid2doid[mim_id]
                for pmid in pmid_list:
                    newrow = [canon, gene_name, hgnc_id, inheritance, disease_name, pmid, mim_id, do_id]
                    row_list.append(newrow)
    return row_list




    return




def extract_glycosylation_motifs_ds(species):
    
    seq_file = path_obj["unreviewed"] +  "%s_protein_canonicalsequences.fasta" % (species)

    newrow = ["uniprotkb_canonical_ac", "start_pos", "end_pos", "motif"] 
    print "\"%s\"" % ("\",\"".join(newrow))
    
    cmd = "perl /software/ps_scan/ps_scan.pl -d %s/prosite/prosite.dat -p PS00001 %s " % (path_obj["downloads"], seq_file)
    
    lines = commands.getoutput(cmd).split("\n")
    for line in lines:
        if line[0:1] == ">":
            canon = line.split("|")[1]
        else:
            parts = line.strip().split(" ")
            newrow = [canon, parts[0], parts[2], parts[4]]
            print "\"%s\"" % ("\",\"".join(newrow))
    return




def extract_isoform_alignments_ds(species):


    work_book = {}
    work_book["masterlist"] = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(work_book["masterlist"], in_file, ",")


    isoform2canon = {}
    for row in work_book["masterlist"]["data"]:
        for isoform in [row[-2], row[-1]]:
            if isoform != "":
                isoform2canon[isoform] = row[0]


    seq_hash = {}
    for fasta_file in glob.glob(path_obj["unreviewed"] +  "%s_protein_allsequences.fasta" % (species)):
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id.split("|")[1]
            seq_hash[seq_id] = str(record.seq.upper())

    seq_set = {}
    for isoform in isoform2canon:
        canon = isoform2canon[isoform]
        aln_file = "alignments/isoformset.uniprotkb.%s.aln" % (canon)
        if os.path.isfile(aln_file) == True:
            file_size = int(os.path.getsize(aln_file))
            if file_size > 0:
                continue
        if canon not in seq_set:
            seq_set[canon] = []
        seq_set[canon].append(">%s\n%s" % (isoform, seq_hash[isoform]))


    for canon in seq_set:
        if len(seq_set[canon]) > 1:
            seq_file = "alignments/isoformset.%s.fasta" % (canon)
            dnd_file = "alignments/isoformset.%s.dnd" % (canon)
            with open(seq_file, "w") as FW:
                FW.write("%s\n" % ("\n\n".join(seq_set[canon])))
            cmd = "/software/clustalw2/clustalw2 -align -type=protein -infile=%s " % (seq_file)
            x = commands.getoutput(cmd)
            cmd = "rm -f %s %s" % (seq_file, dnd_file)
            x = commands.getoutput(cmd)
            print "aligned isoformset for canon=%s" % (canon)


    return




def extract_pubchem_cid_smiles_canonical_sequences_ds(species):


    work_book = {}
    work_book["masterlist"] = {}
    idmap_file = path_obj["unreviewed"] +  "%s_glycan_masterlist.csv" % (species)
    libgly.load_sheet(work_book["masterlist"], idmap_file, ",")

    cid2glytoucan = {}
    for row in work_book["masterlist"]["data"]:
        if row[3][0:3] == "CID":
            cid = row[3][3:]
            cid2glytoucan[cid] = row[0]

    row = ["glytoucan_ac","pubchem_cid","smiles_canonical"]
    print "\"%s\"" % ("\",\"".join(row))
    in_file = path_obj["downloads"] + "pubchem/compound/cid2smiles.tsv"
    with open(in_file, "r") as FR:
        for line in FR:
            cid, smiles = line.strip().split()
            if cid in cid2glytoucan:
                row = [cid2glytoucan[cid], cid, smiles]
                print "\"%s\"" % ("\",\"".join(row))



    return



def extract_expression_disease_ds(species, in_file):

    work_book = {}
    sheet_name = "masterlist"
    work_book[sheet_name] = {}
    idmap_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")


    ac2canon = {}
    for row in work_book[sheet_name]["data"]:
        canon = row[0]
        ac = canon.split("-")[0]
        ac2canon[ac] = canon


    sheet_name = "diseaseexpression"
    work_book[sheet_name] = {}
    libgly.load_sheet(work_book[sheet_name], in_file, ",")


    row = ["uniprotkb_canonical_ac","significance","direction", "do_id", "do_name", "parent_doid", "parent_doname"]
    print "\"%s\"" % ("\",\"".join(row))

    f_list = work_book[sheet_name]["fields"]
    for row in work_book[sheet_name]["data"]:
        if row[0] not in ac2canon:
            continue
        canon = ac2canon[row[0]]

        if row[f_list.index("doid")] in ["3963", "0070003", "3119"]:
            continue
        parent_doid, parent_doname = "", ""
        if len(row[f_list.index("parent_doname")]) > 0:
            parent_doid, parent_doname = row[f_list.index("parent_doid")], row[f_list.index("parent_doname")]
        newrow = [
            canon,
            row[f_list.index("significance")],
            row[f_list.index("direction")],
            row[f_list.index("doid")],
            row[f_list.index("doname")],
            parent_doid, 
            parent_doname
        ]
        print "\"%s\"" % ("\",\"".join(newrow))




def extract_expression_normal_ds(species, in_file):

    work_book = {}
    sheet_name = "masterlist"
    work_book[sheet_name] = {}
    idmap_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")
   

    ac2canon = {}
    for row in work_book[sheet_name]["data"]:
        canon = row[0]
        ac = canon.split("-")[0]
        ac2canon[ac] = canon

    sheet_name = "normalexpression"
    work_book[sheet_name] = {}
    libgly.load_sheet(work_book[sheet_name], in_file, ",")        

    row = ["uniprotkb_canonical_ac","score","uberon_dev_id","sex", "uberon_anatomy_id","expression_call","uberon_name"]
    print "\"%s\"" % ("\",\"".join(row))
    f_list = work_book[sheet_name]["fields"]
    for row in work_book[sheet_name]["data"]:
        if row[0] not in ac2canon:
            continue
        canon = ac2canon[row[0]]
        newrow = [
            canon, 
            row[f_list.index("expressionScore")], 
            row[f_list.index("uberonDevelopmentId")],
            row[f_list.index("sex")],
            row[f_list.index("uberonAnatomyId")],
            row[f_list.index("expressionCall")],
            row[f_list.index("uberon_name")]
        ]
        print "\"%s\"" % ("\",\"".join(newrow))






def extract_mutation_ds(species):



    work_book = {}
    sheet_name = "masterlist"
    work_book[sheet_name] = {}
    idmap_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")

    canon_list = []
    for row in work_book[sheet_name]["data"]:
        canon_list.append(row[0])

       
    row = ["uniprotkb_canonical_ac","aa_pos","ref_aa","alt_aa","chromosome_id","chr_pos","ref_nt","alt_nt", "patients_positive", "patients_tested", "mut_freq","data_source","do_id","do_name"]
    print "\"%s\"" % ("\",\"".join(row))

    
    in_file = path_obj["downloads"] + "biomuta/biomuta.csv"
    f_list = []
    freq_dist = {}
    cutoff_one = 0.01
    cutoff_two = 0 
    with open(in_file, "r") as FR:
        csv_grid = csv.reader(FR, delimiter=",", quotechar='\"')
        row_count = 0
        for row in csv_grid:
            row_count += 1
            if row_count == 1:
                f_list = row
                continue
            isoform = row[f_list.index("uniprot_canonical_ac")]
            patients_tested = int(row[f_list.index("patients_tested")])
            patients_positive = int(row[f_list.index("patients_positive")])
            mut_freq = float(row[f_list.index("mut_freq")])
            if patients_positive < cutoff_two or float(mut_freq) < cutoff_one:
                continue
            if isoform not in canon_list:
                continue
            if row[f_list.index("do_id")] in ["3963", "0070003", "3119"]:
                continue
            if row[f_list.index("ref_aa")] == row[f_list.index("alt_aa")]:
                continue
            newrow = [
                isoform,
                row[f_list.index("aa_pos")],
                row[f_list.index("ref_aa")],
                row[f_list.index("alt_aa")],
                row[f_list.index("chr_id")],
                row[f_list.index("chr_pos")],
                row[f_list.index("ref_nt")],
                row[f_list.index("alt_nt")],
                row[f_list.index("patients_positive")],
                row[f_list.index("patients_tested")],
                row[f_list.index("mut_freq")],
                row[f_list.index("source")],
                row[f_list.index("do_id")],
                row[f_list.index("do_name")]
            ]
            print "\"%s\"" % ("\",\"".join(newrow))



    return


def run_stretcher(seq_hash, seqid_one, seqid_two):

    cmd = "rm /tmp/seq_one.fasta /tmp/seq_two.fasta /tmp/stretcher.aln"
    x = commands.getoutput(cmd)

    with open("/tmp/seq_one.fasta", "w") as FW:
        FW.write("%s\n" % (">" + seqid_one + "\n" + seq_hash[seqid_one]))

    with open("/tmp/seq_two.fasta", "w") as FW:
        FW.write("%s\n" % (">" + seqid_two + "\n" + seq_hash[seqid_two]))

    cmd = "stretcher %s %s %s" % ("/tmp/seq_one.fasta", "/tmp/seq_two.fasta", "/tmp/stretcher.aln")
    x = commands.getoutput(cmd)
    
    cmd = 'grep "# Identity:" /tmp/stretcher.aln'
    x = commands.getoutput(cmd).strip().split("(")[-1].split("%")[0]
    return x



def extract_genelocus_ds(species):

    gtf_file = path_obj["downloads"] + "ucsc/gtf/%s.gtf" % (species)

    gene_locus = {}
    with open(gtf_file, 'r') as FR:
        data_frame = csv.reader(FR, delimiter='\t', quotechar='|')
        rowCount = 0
        for row in data_frame:
            if row[0][0] == "#":
                continue
            rowCount += 1
            if row[2] == "gene":
                att_list = row[-1].split(";")
                att_dict = {}
                for att in att_list:
                    if att.strip() == "":
                        continue
                    name, value = att.strip().split(" ")
                    att_dict[name] = value.replace("\"", "")
                if "gene_name" in att_dict:
                    strand = ""
                    strand = "0" if row[6].strip() == "-" else strand
                    strand = "1" if row[6].strip() == "+" else strand
                    newrow = [att_dict["gene_id"],row[0],row[3],row[4],strand]
                    gene_locus[att_dict["gene_name"]] = newrow


    sheet_obj = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(sheet_obj, in_file, ",")
    f_list = sheet_obj["fields"]
    newrow = ["uniprotkb_canonical_ac","gene_symbol", "ensembl_gene_id", "chromosome_id", 
            "start_pos", "end_pos", "strand"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    for row in sheet_obj["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        if canon not in seen:
            gene_name = row[f_list.index("gene_name")]
            if gene_name in gene_locus:
                newrow = [canon, gene_name] + gene_locus[gene_name]
                print "\"%s\"" % ("\",\"".join(newrow))
        seen[canon] = True

    
    return


def extract_binary_interactions_ds(species):
    
    data_grid = { "ac2canon":{}, "isoformlist":{}, "interaction":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_interaction(data_grid, species)

    newrow = ["uniprotkb_canonical_ac","intact_ac","participant_uniprotkb_ac","participant_intact_ac",
            "participant_uniprotkb_id","participant_gene_symbol", "participant_taxid","experiments"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for ac in data_grid["interaction"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            for o in data_grid["interaction"][ac]:
                newrow = [canon,o["intactac"], o["puniprotac"],o["pintactac"],o["puniprotid"], 
                        o["pgenename"], o["ptaxid"],o["experiments"]]
                print "\"%s\"" % ("\",\"".join(newrow))
    return



def extract_enzyme_annotation_uniprotkb_ds(species):

    data_grid = { "ac2canon":{}, "isoformlist":{}, "enzymeinfo":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_enzymelist_one(data_grid, species)
    sparqlutil.load_enzymelist_two(data_grid, species)

    newrow = ["uniprotkb_canonical_ac", "enzyme_ec", "enzyme_activity"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for ac in data_grid["enzymeinfo"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            for o in data_grid["enzymeinfo"][ac]:
                newrow = [canon,o["ec"],o["activity"]]
                print "\"%s\"" % ("\",\"".join(newrow))



    return




def extract_glycoenzymes_ds(species, dataset):

    data_grid = { "ac2canon":{}, "isoformlist":{}}
    sparqlutil.load_isoformlist(data_grid, species)



    work_book = {}
    sheet_name = "masterlist"
    work_book[sheet_name] = {}
    idmap_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")

    seen_canon = {}
    for row in work_book[sheet_name]["data"]:
        canon = row[0]
        ac = canon.split("-")[0]
        seen_canon[canon] = True


    sheet_name = dataset
    work_book[sheet_name] = {}
    in_file = "compiled/%s_protein_%s.csv" % (species, dataset)
    libgly.load_sheet(work_book[sheet_name], in_file, ",")

    f_list = work_book[sheet_name]["fields"]
    f_list[0] = "uniprotkb_canonical_ac"
    print "\"%s\"" % ("\",\"".join(f_list))
    ncols = len(f_list)

    for row in work_book[sheet_name]["data"]:
        ac = row[0]
        if ac not in data_grid["ac2canon"]:
            continue
        canon = data_grid["ac2canon"][ac]
        if canon not in seen_canon:
            continue
        if len(row) != ncols:
            continue
        row[0] = canon
        print "\"%s\"" % ("\",\"".join(row))

    return





def extract_function_uniprotkb_ds(species):

    data_grid = { "ac2canon":{}, "isoformlist":{},  "function":{}}

    sparqlutil.load_function(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)
    
    row = ["uniprotkb_canonical_ac","database_name", "database_id", "evidence","annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    for ac in data_grid["function"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            ann_list = data_grid["function"][ac] if ac in data_grid["function"] else []
            for ann in ann_list:
                word_list = ann.split(" ")
                pmid_list = []
                for word in word_list:
                    if word.find("PubMed:") != -1:
                        word = word.split(")")[0]
                        word = word.replace("(", "").replace(")", "").replace(",", "")
                        pmid_list.append(word[7:])
                if pmid_list == []:
                    pmid_list = [""]
                for pmid in sorted(set(pmid_list)):
                    ann = ann.replace("\"", "`")
                    row = [canon, "UniProtKB", ac, pmid, ann]
                    print "\"%s\""  % ("\",\"".join(row))



    
    return


    





def extract_disease_ds(species):

    data_grid = {
        "ac2canon":{}, 
        "isoformlist":{}, 
        "ac2omim":{}, 
        "ac2mondo":{}
    }

    glycophenotype_sheet = {}
    in_file = path_obj["downloads"] + "ohsu/human_glycophenotype.csv"
    libgly.load_sheet(glycophenotype_sheet,in_file, ",")

    f_list = glycophenotype_sheet["fields"]
    for row in glycophenotype_sheet["data"]:
        ac = row[f_list.index("uniprotkb_ac")]
        mondo_id = row[f_list.index("mondo_id")].split(":")[1]
        mondo_label = row[f_list.index("mondo_label")]
        if ac not in data_grid["ac2mondo"]:
            data_grid["ac2mondo"][ac] = []
        data_grid["ac2mondo"][ac].append({"id":mondo_id, "label":mondo_label})

    sparqlutil.load_disease_mim(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)


    seen = {}
    row = ["uniprotkb_canonical_ac","database_name", "database_id", "database_label"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        if ac in data_grid["ac2mondo"]:
            for o in data_grid["ac2mondo"][ac]:
                row = [canon, "mondo", o["id"], o["label"]]
                row_str = ",".join(row)
                if row_str not in seen:
                    print "\"%s\""  % ("\",\"".join(row))
                    seen[row_str] = True
        if ac in data_grid["ac2omim"]:
            for mim_id in data_grid["ac2omim"][ac]:
                row = [canon, "omim", mim_id, ""]
                row_str = ",".join(row)
                if row_str not in seen:
                    print "\"%s\""  % ("\",\"".join(row))
                    seen[row_str] = True

    row_list = get_genomics_england_disease_rowlist(species)
    for row in row_list[1:]:
        #make row from gene_name
        newrow = [row[0], "genomics_england", row[1],row[4]]
        newrow_str = ",".join(newrow)
        if newrow_str not in seen:
            print "\"%s\""  % ("\",\"".join(newrow))
            seen[newrow_str] = True
        
        #make row from mim_id
        newrow = [row[0], "omim", row[6],""]
        newrow_str = ",".join(newrow)
        if newrow_str not in seen:
            print "\"%s\""  % ("\",\"".join(newrow))
            seen[newrow_str] = True

 
        mim_id




    return


def extract_citations_reactome_ds(species):

    black_list = get_blacklisted_pmids(species)

    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_reactome.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    pathway2canon = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        pathway_id = row[f_list.index("database_id")]
        if pathway_id not in pathway2canon:
            pathway2canon[pathway_id] = []
        if canon not in pathway2canon[pathway_id]:
            pathway2canon[pathway_id].append(canon)

    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_reactions_reactome.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    newrow = ["uniprotkb_canonical_ac","pmid","title","journal_name","publication_date", "authors"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    for row in data_frame["data"]:
        pmid = row[f_list.index("pmid")]
        pathway_id = row[f_list.index("pathway_id")]
        if pmid in black_list:
            continue
        if pathway_id not in pathway2canon:
            continue
        newrow = get_citation(pmid)
        if newrow != []:
            for canon in pathway2canon[pathway_id]:
                combo_id = "%s %s" % (canon, pmid)
                if combo_id not in seen:
                    print "\"%s\"" % ("\",\"".join([canon] + newrow))
                    seen[combo_id] = True



    return


def get_blacklisted_pmids(species):

    black_list = []
    in_file = "compiled/%s_protein_blacklisted_pmids_uniprotkb.csv" % (species)
    if os.path.isfile(in_file) == True:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        for row in data_frame["data"]:
            black_list.append(row[0])
        black_list = sorted(set(black_list))
       
    return black_list



def extract_citations_uniprotkb_ds(species):

    black_list = get_blacklisted_pmids(species)

    data_grid = {"ac2canon":{}, "isoformlist":{}, "citelist":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_citelist(data_grid, species)

    row = ["uniprotkb_canonical_ac","pmid","title","journal_name", "publication_date", "authors"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["citelist"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            for obj in data_grid["citelist"][ac]:
                if obj["pmid"] in black_list:
                    continue
                authors = ", ".join(obj["authorlist"])
                row = [canon,obj["pmid"], obj["journaltitle"], obj["journalname"], obj["pubdate"], authors]
                print "\"%s\""  % ("\",\"".join(row))



def extract_canonicalsequences_ds(species):

    data_grid = {"ac2canon":{}, "genename":{}, "recnames":{}, "proteinid":{}, "isoformlist":{}, "isoforminfo":{}, "isoformseq":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)
    sparqlutil.load_isoformseq(data_grid, species)
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_recnames(data_grid, species)
    sparqlutil.load_genename(data_grid, species)

    tax_name = species_obj[species]["taxname"]


    for ac in data_grid["isoformlist"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            protein_id = data_grid["proteinid"][ac] if ac in data_grid["proteinid"] else ""
            gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
            full_name = data_grid["recnames"][ac]["fullname"] if ac in data_grid["recnames"] else ""
            
            seq = data_grid["isoformseq"][canon]
            status = data_grid["isoforminfo"][canon]["reviewed"]
            ac_lbl = "tr|%s|%s" % (canon,protein_id)
            desc = "%s OS=%s GN=%s" %(full_name, tax_name, gene_name)
            if status == "1":
                ac_lbl = "sp|%s|%s" % (canon,protein_id)
            seq_obj = SeqRecord(Seq(seq,IUPAC.protein),id=ac_lbl, name=ac_lbl, description=desc)
            print "%s" % (seq_obj.format("fasta"))




def extract_allsequences_ds(species):
        
    data_grid = {"ac2canon":{}, "genename":{}, "recnames":{}, "proteinid":{}, "isoformlist":{}, "isoforminfo":{}, "isoformseq":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)
    sparqlutil.load_isoformseq(data_grid, species)
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_recnames(data_grid, species)
    sparqlutil.load_genename(data_grid, species)

    tax_name = species_obj[species]["taxname"]

    for ac in data_grid["isoformlist"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            protein_id = data_grid["proteinid"][ac] if ac in data_grid["proteinid"] else ""
            gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
            full_name = data_grid["recnames"][ac]["fullname"] if ac in data_grid["recnames"] else ""
            for isoform in data_grid["isoformlist"][ac]:
                seq = data_grid["isoformseq"][isoform]
                status = data_grid["isoforminfo"][isoform]["reviewed"] 
                ac_lbl = "tr|%s|%s|%s" % (isoform,canon,protein_id)
                desc = "%s OS=%s GN=%s" %(full_name, tax_name, gene_name)
                if status == "1":
                    ac_lbl = "sp|%s|%s|%s" % (isoform,canon,protein_id)
                seq_obj = SeqRecord(Seq(seq,IUPAC.protein),id=ac_lbl, name=ac_lbl, description=desc)
                print "%s" % (seq_obj.format("fasta"))
                                            



    return

def extract_info_uniprotkb_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "proteinid":{}, "isoformmass":{}, "isoformlen":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_isoformmass(data_grid, species)
    sparqlutil.load_isoformlen(data_grid, species)


    row = ["uniprotkb_canonical_ac","uniprotkb_id","uniprotkb_protein_mass","uniprotkb_protein_length"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        protein_id = data_grid["proteinid"][ac] if ac in data_grid["proteinid"] else ""
        canon_mass = data_grid["isoformmass"][canon] if canon in data_grid["isoformmass"] else -1
        canon_len = data_grid["isoformlen"][canon] if canon in data_grid["isoformlen"] else -1

        row = [canon, protein_id, canon_mass, canon_len]
        print "\"%s\""  % ("\",\"".join(row))




    return

def extract_ac2pdb_ds(species):

    ds_name = "ac2pdb"

    data_grid = {"ac2canon":{}, "isoformlist":{}, "ac2pdb":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_ac2pdb(data_grid, species)
    
    row = ["uniprotkb_canonical_ac","pdb_id"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        if ac in data_grid[ds_name]:
            for pdb_id in data_grid[ds_name][ac]:
                row = [canon, pdb_id]
                print "\"%s\""  % ("\",\"".join(row))


def extract_xrefs_ds(species, ds_name):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "reactome":{}}
    sparqlutil.load_isoformlist(data_grid, species)

    sheet_obj = {}
    sparqlutil.load_ac2xref(sheet_obj, species, config_obj["xref"][ds_name])

    row = ["uniprotkb_canonical_ac","database_id","database_label"]
    print "\"%s\""  % ("\",\"".join(row))
    seen = {}
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        if ac in sheet_obj:
            #consider one-to-one mapping only for xref_refseq
            if ds_name == "xref_refseq" and len(sheet_obj[ac]) > 1:
                continue
            for o in sheet_obj[ac]:
                row = [canon, o["id"], o["label"].encode('ascii', 'ignore').decode('ascii')]
                print "\"%s\""  % ("\",\"".join(row))
                seen[canon] = True


    #Add isoform level mappings for xref_refseq
    if ds_name == "xref_refseq":
        sheet_obj = {}
        sparqlutil.load_refseq2isoform(sheet_obj, species)
        for ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            if canon in seen:
                continue
            if canon in sheet_obj:
                selected = sheet_obj[canon][0]
                for refseq_isoform in sheet_obj[canon]:
                    if refseq_isoform[0:3] == "NP_":
                        selected = refseq_isoform
                        break
                row = [canon, selected, ""]
                print "\"%s\""  % ("\",\"".join(row))




    return



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

def extract_citations_refseq_ds(species):

    black_list = get_blacklisted_pmids(species)

    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_function_refseq.csv" % (species)
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



def extract_recnames_ds(species):

    ds_name = "recnames"

    data_grid = {"ac2canon":{}, "isoformlist":{}, "recnames":{}, "shortname":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_recnames(data_grid, species)

    row = ["uniprotkb_canonical_ac","recommended_name_full","recommended_name_short"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        full_name, short_name = "", ""
        if ac in data_grid[ds_name]:
            obj = data_grid[ds_name][ac]
            full_name = obj["fullname"] if "fullname" in obj else ""
            short_name = obj["shortname"] if "shortname" in obj else ""
            row = [canon, full_name, short_name]
            print "\"%s\""  % ("\",\"".join(row))



def extract_altnames_ds(species):

    ds_name = "altnames"
    data_grid = {"ac2canon":{}, "isoformlist":{}, "altnames":{}, "shortname":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_altnames(data_grid, species)

    row = ["uniprotkb_canonical_ac","alternative_name_full","alternative_name_short"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        full_name, short_name = "", ""
        if ac in data_grid["altnames"]:
            for name_uri in data_grid["altnames"][ac]:
                full_name, short_name = "", ""
                if "fullname" in data_grid["altnames"][ac][name_uri]:
                    full_name = data_grid["altnames"][ac][name_uri]["fullname"]
                if "shortname" in data_grid["altnames"][ac][name_uri]:
                    short_name = data_grid["altnames"][ac][name_uri]["shortname"] 
                row = [canon, full_name, short_name]
                print "\"%s\""  % ("\",\"".join(row))



    return


def extract_pdb_shortlist_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "pdbinfo":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_pdbinfo(data_grid, species)
    row = ["uniprotkb_canonical_ac","pdb_id", "method", "resolution"]
    print "\"%s\""  % ("\",\"".join(row))

    seen = {}
    selected_dict = {"xray":{}, "nmr":{}}
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        if ac in data_grid["pdbinfo"]:
            for obj in data_grid["pdbinfo"][ac]:
                if canon not in selected_dict["xray"]:
                    selected_dict["xray"][canon] = []
                if canon not in selected_dict["nmr"]:
                    selected_dict["nmr"][canon] = []
                if obj["method"] == "X-Ray_Crystallography":
                    selected_dict["xray"][canon].append({"pdbid":obj["pdbid"], "res":obj["res"]})
                elif obj["method"] == "Electron_Microscopy":
                    selected_dict["nmr"][canon].append({"pdbid":obj["pdbid"], "res":obj["res"]})

    seen = {}
    #select pdbid with X-ray method and minimum resolution
    for canon in selected_dict["xray"]:
        min_res = 1000.0
        pdb_id = ""
        for obj in selected_dict["xray"][canon]:
            pdb_id = obj["pdbid"] if obj["res"] < min_res else pdb_id
            min_res = obj["res"] if obj["res"] < min_res else min_res
        if pdb_id == "":
            continue
        row = [canon,pdb_id,"X-Ray_Crystallography", str(min_res)]
        row_str = "\",\"".join(row)
        print "\"%s\""  % (row_str)
        seen[canon] = True
    
    for canon in selected_dict["nmr"]:
        if canon in seen:
            continue
        obj = selected_dict["nmr"][canon][0]
        pdb_id = obj["pdbid"]
        min_res = obj["res"]
        #for obj in selected_dict["nmr"][canon]:
        #    print canon, obj
        if pdb_id == "":
            continue
        row = [canon,pdb_id, "Electron_Microscopy", str(min_res)]
        row_str = "\",\"".join(row)
        print "\"%s\""  % (row_str)
        seen[canon] = True



    return



def extract_ptm_annotation_uniprotkb_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "ptmann":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_ptm_annotation(data_grid, species)
    row = ["uniprotkb_canonical_ac","pmid", "eco_id", "ptm_annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    seen = {}
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        if ac in data_grid["ptmann"]:
            for obj in data_grid["ptmann"][ac]:
                row = [canon,obj["pmid"], obj["ecoid"], obj["ann"]]
                row_str = "\",\"".join(row)
                if row_str not in seen:
                    print "\"%s\""  % (row_str)
                    seen[row_str] = True

    return


def extract_ntdata_ds(species):
    sp = species_obj[species]["taxname"].lower().replace(" ", "-")
    cmd = "cp downloads/ebi/current/uniprot-proteome-%s.nt unreviewed/%s_protein_ntdata.nt" % (sp, species)
    x = commands.getoutput(cmd)

    return



def extract_site_annotation_uniprotkb_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "genename":{}, "siteann":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_genename(data_grid, species)

    ann_type_list = [
        "Active_Site_Annotation", 
        "Binding_Site_Annotation", 
        "Site_Annotation",
        "Calcium_Binding_Annotation",
        "Disulfide_Bond_Annotation",
        "Glycosylation_Annotation",
        "Metal_Binding_Annotation",
        "Modified_Residue_Annotation",
        "NP_Binding_Annotation",
        "Nucleotide_Binding_Annotation",
        "PTM_Annotation",
        "Signal_Peptide_Annotation"
    ]
    for ann_type in ann_type_list:
        sparqlutil.load_site_annotation(data_grid, species, ann_type)
   
    row = ["uniprotkb_canonical_ac","gene_symbol", "ann_type", "start_pos","end_pos", "annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        if ac in data_grid["siteann"]:
            for obj in data_grid["siteann"][ac]:
                row = [canon, gene_name,obj["anntype"], str(obj["start"]), str(obj["end"]),obj["ann"]]
                print "\"%s\""  % ("\",\"".join(row))
    return

def extract_go_annotation_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "genename":{}, "goann":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_go_annotation(data_grid, species)
    sparqlutil.load_genename(data_grid, species)


    row = ["uniprotkb_canonical_ac","gene_symbol", "go_term_id","go_term_label", "go_term_category"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        if ac in data_grid["goann"]:
            for obj in data_grid["goann"][ac]:
                row = [canon, gene_name, obj["goid"], obj["goterm"],obj["gocat"]]
                print "\"%s\""  % ("\",\"".join(row))

    return




def extract_transcriptlocus_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "isoforminfo":{}, "locusinfo":{}}
    sparqlutil.load_locusinfo(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)

    row = ["uniprotkb_canonical_ac","uniprotkb_isoform_ac","transcript_id","peptide_id","chromosome_id",
            "start_pos","end_pos","strand"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["isoformlist"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            for isoform in data_grid["isoformlist"][ac]:
                if isoform in data_grid["locusinfo"]:
                    o = data_grid["locusinfo"][isoform]
                    row = [canon,isoform,o["trsid"], o["pepid"],
                            str(o["chrid"]),str(o["startpos"]),str(o["endpos"]),str(o["strand"])]
                    print "\"%s\"" % ("\",\"".join(row))




def extract_masterlist_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "genename":{}, "isoforminfo":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_genename(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)



    row = ["uniprotkb_canonical_ac","status","gene_name","reviewed_isoforms","unreviewed_isoforms"]
    print "\"%s\""  % ("\",\"".join(row))

    for ac in data_grid["isoformlist"]:
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        list_one, list_two = [], []
        for isoform in data_grid["isoformlist"][ac]:
            if data_grid["isoforminfo"][isoform]["reviewed"] == "1":
                list_one.append(isoform)
            if data_grid["isoforminfo"][isoform]["reviewed"] == "0":
                list_two.append(isoform)
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            canon_status = "reviewed" if data_grid["isoforminfo"][canon]["reviewed"] == "1" else "unreviewed"
            n = len(list_one) if len(list_one) > len(list_two) else len(list_two)
            for j in xrange(0, n):
                val_one = list_one[j] if j < len(list_one) else ""
                val_two = list_two[j] if j < len(list_two) else ""
                row = [canon, canon_status, gene_name, val_one, val_two]
                print "\"%s\""  % ("\",\"".join(row))

    return

def extract_info_refseq_ds(species, tax_id):

    refseq2canon = {}
    canon2refseq = {}
    map_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
    with open(map_file, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_frame:
            row_count += 1
            canon, refseq = row[0], row[1]
            if canon not in canon2refseq:
                canon2refseq[canon] = refseq
                refseq2canon[refseq] = canon
 

    data_frame = {}
    field = "xxx"
    in_file = path_obj["downloads"] + "ncbi/refseq/refseq_protein_all_%s.gpff" % (tax_id)
    with open(in_file, "r") as FR:
        for line in FR:
            newfield = line[0:12].strip() 
            field = newfield if len(newfield) > 0 else field
            value = line[12:].strip()
            if field == "VERSION":
                ac = line.split(" ")[-1].strip()
                data_frame[ac] = {"ac":ac, "summary":""}
                flag = False
            elif field == "COMMENT":
                if line.strip()[0:8] == "Summary:":
                    flag = True
                if line.strip() == "":
                    flag = False
                if flag:
                    data_frame[ac]["summary"] += line.strip() + " "


    row = ["uniprotkb_canonical_ac","p_refseq_ac_best_match","refseq_protein_name", "refseq_protein_length","refseq_protein_summary"]
    print "\"%s\"" % ("\",\"".join(row))
    for rec in SeqIO.parse(in_file, "genbank"):
        ac = rec.id
        for feat in rec.features:
            if feat.type in ["Protein"]:
                summary = data_frame[ac]["summary"]
                summary = summary.replace("\"", "`")
                summary = summary.replace("Summary: ", "")
                product = feat.qualifiers["product"][0]
                product = product.replace("\"", "`")
                seq_len = len(str(rec.seq))
                if ac in refseq2canon:
                    uniprotkb_canonical_ac = refseq2canon[ac]
                    row = [uniprotkb_canonical_ac, ac, product, str(seq_len), summary]
                    print "\"%s\"" % ("\",\"".join(row))



def extract_function_refseq_ds(species, tax_id):


    refseq2canon = {}
    canon2refseq = {}
    map_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
    with open(map_file, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_frame:
            row_count += 1
            canon, refseq = row[0], row[1]
            if canon not in canon2refseq:
                canon2refseq[canon] = refseq
                refseq2canon[refseq] = canon
   

    in_file = path_obj["downloads"] + "ncbi/refseq/refseq_protein_all_%s.gpff" % (tax_id)
    data_frame = {}
    field = "xxx"
    with open(in_file, "r") as FR:
        for line in FR:
            newfield = line[0:12].strip() 
            field = newfield if len(newfield) > 0 else field
            value = line[12:].strip()
            if field == "VERSION":
                ac = line.split(" ")[-1].strip()
                data_frame[ac] = {"ac":ac, "references":[]}
            elif field == "REFERENCE":
                data_frame[ac]["references"].append({"pubmed":"", "remark":""})
            elif field == "REMARK":
                data_frame[ac]["references"][-1]["remark"] += value.strip() + " "
            elif field == "PUBMED":
                data_frame[ac]["references"][-1]["pubmed"] = value


    row = ["uniprotkb_canonical_ac","database_name", "database_id", "evidence","annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    for ac in data_frame:
        for o in data_frame[ac]["references"]:
            if o["remark"].find("GeneRIF:") != -1:
                o["remark"] = o["remark"][8:].strip()
                o["remark"] = o["remark"][0].upper() + o["remark"][1:]
                o["remark"] = o["remark"].replace("\"", "`")
                uniprotkb_canonical_ac = "N/A"
                if ac in refseq2canon:
                    uniprotkb_canonical_ac = refseq2canon[ac]
                    row = [uniprotkb_canonical_ac, "RefSeq", ac, o["pubmed"], o["remark"]]
                    print "\"%s\"" % ("\",\"".join(row))

    
    return




def extract_participants_reactome_ds(species):

    sheet_obj = {"participants":{}}
    sparqlutil.load_participants_reactome(sheet_obj, [species])

    tmprow = ["reaction_id","participant_id", "participant_name","role", "xref_id", "xref_type"]
    print "\"%s\""  % ("\",\"".join(tmprow))
    
    for reaction_id in sheet_obj["participants"]:
        for role in ["input", "output"]:
            for o in sheet_obj["participants"][reaction_id][role]:
                name = o["name"].encode('ascii', 'ignore').decode('ascii')
                tmprow = [
                    reaction_id,o["id"], name, role, o["xrefid"], o["xreftype"]
                ]
                print "\"%s\""  % ("\",\"".join(tmprow))

    return



def extract_reactions_reactome_ds(species):

    sheet_obj = {"reactions":[]}
    sparqlutil.load_reactions_reactome(sheet_obj, [species])

    tmprow = ["reaction_id","reaction_name","cellular_location", "pmid", "pathway_id"]
    print "\"%s\""  % ("\",\"".join(tmprow))

    for o in sheet_obj["reactions"]:
        tmprow = [
            o["reactionid"],
            o["reactionname"].encode('ascii', 'ignore').decode('ascii'),
            o["cellularlocation"].encode('ascii', 'ignore').decode('ascii'),
            o["pmid"],
            o["pathwayid"]
        ]
        print "\"%s\""  % ("\",\"".join(tmprow))

    return


def extract_misc_ds(dataset, species):

    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        isoform_list = [row[f_list.index("reviewed_isoforms")], row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon
            ac2canon[isoform] = canon
            ac2canon[canon] = canon


    if dataset in ["cleaved_signal_peptide_sequence", "signal_peptide_sequence"]:
        in_file = "compiled/%s_protein_%s.fasta" % (species, dataset)
        seq_hash = {}
        log_file_one = path_obj["logs"] +  "%s_protein_%s.1.log" % (species, dataset)
        FL1 = open(log_file_one, "w")
        for record in SeqIO.parse(in_file, "fasta"):
            canon = record.id.split("|")[1]
            seq = str(record.seq.upper())
            if canon not in ac2canon:
                FL1.write(">%s\n%s\n" % (canon, seq) )
            else:
                print ">%s\n%s" % (canon, seq)
        FL1.close()
    else:
        data_frame = {}
        in_file = "compiled/%s_protein_%s.csv" % (species, dataset)
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        print "\"%s\""  % ("\",\"".join(["uniprotkb_canonical_ac"] + f_list[1:]))
    
        log_file_one = path_obj["logs"] +   "%s_protein_%s.1.log" % (species, dataset)
        FL1 = open(log_file_one, "w")
        for row in data_frame["data"]:
            ac = row[0]
            if ac in ac2canon:
                canon = ac2canon[ac]
                newrow = [canon] + row[1:]
                print "\"%s\""  % ("\",\"".join(newrow))
            else:
                FL1.write("\"%s\"\n" % ("\",\"".join(row)))
        FL1.close()

    return


###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")
    parser.add_option("-d","--dataset",action="store",dest="dataset",help="[masterlist, transcriptlocus]")



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
    if dataset == "masterlist":
        extract_masterlist_ds(species)
    elif dataset == "transcriptlocus":
        extract_transcriptlocus_ds(species)
    elif dataset == "recnames":
        extract_recnames_ds(species)
    elif dataset == "altnames":
        extract_altnames_ds(species)
    elif dataset == "ac2pdb":
        extract_ac2pdb_ds(species)
    elif dataset == "info_uniprotkb":
        extract_info_uniprotkb_ds(species)
    elif dataset == "allsequences":
        extract_allsequences_ds(species)
    elif dataset == "canonicalsequences":
        extract_canonicalsequences_ds(species)
    elif dataset == "citations_uniprotkb":
        extract_citations_uniprotkb_ds(species)
    elif dataset == "citations_reactome":
        extract_citations_reactome_ds(species)
    elif dataset == "disease":
        extract_disease_ds(species)
    elif dataset == "function_uniprotkb":
        extract_function_uniprotkb_ds(species)
    elif dataset in config_obj["xref"]:
        extract_xrefs_ds(species, dataset)
    elif dataset == "go_annotation":
        extract_go_annotation_ds(species)
    elif dataset == "site_annotation_uniprotkb":
        extract_site_annotation_uniprotkb_ds(species)
    elif dataset == "ptm_annotation_uniprotkb":
        extract_ptm_annotation_uniprotkb_ds(species)
    elif dataset == "pdb_shortlist":
        extract_pdb_shortlist_ds(species)
    elif dataset == "mutation":
        extract_mutation_ds(species)
    elif dataset == "expression_disease":
        in_file = path_obj["downloads"] + "bioxpress/bioxpress_disease.csv"
        extract_expression_disease_ds(species, in_file)
    elif dataset == "expression_normal":
        in_file = path_obj["downloads"] + "bioxpress/bioxpress_normal.csv"
        extract_expression_normal_ds(species, in_file)
    elif dataset in ["glycohydrolase", "glycosyltransferase"]:
        extract_glycoenzymes_ds(species, dataset)
    elif dataset == "genelocus":
        extract_genelocus_ds(species)
    elif dataset == "pubchem_cid_smiles_canonical_sequences":
        extract_pubchem_cid_smiles_canonical_sequences_ds(species)
    elif dataset == "isoform_alignments":
        extract_isoform_alignments_ds(species)
    elif dataset == "info_refseq":
        extract_info_refseq_ds(species, species_obj[species]["taxid"])
    elif dataset == "function_refseq":
        extract_function_refseq_ds(species, species_obj[species]["taxid"])
    elif dataset == "citations_refseq":
        extract_citations_refseq_ds(species)
    elif dataset == "reactions_reactome":
        extract_reactions_reactome_ds(species)
    elif dataset == "participants_reactome":
        extract_participants_reactome_ds(species)
    elif dataset == "glycosylation_motifs":
        extract_glycosylation_motifs_ds(species)
    elif dataset in ["cdg_disease", "glyco_disease", "proteinwoutsignalp", 
            "cleaved_signal_peptide_sequence", "signal_peptide_sequence"]:
        extract_misc_ds(dataset, species)
    elif dataset == "enzyme_annotation_uniprotkb":
        extract_enzyme_annotation_uniprotkb_ds(species)
    elif dataset == "binary_interactions":
        extract_binary_interactions_ds(species)
    elif dataset == "genomics_england_disease":
        extract_genomics_england_disease_ds(species)
    elif dataset == "ntdata":
        extract_ntdata_ds(species)
















if __name__ == '__main__':
	main()
