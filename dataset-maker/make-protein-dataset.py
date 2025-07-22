import os,sys
import json
import csv
import time

#import requests

from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import GenBank


import commands
import glob
import re
import sparqlutil


import libgly


__version__="1.0"
__status__ = "Dev"


def cluster_pdb_ids(in_dict):

    obj_dict = {}
    for pdb_id in in_dict:
        for r in in_dict[pdb_id]:
            s, e = int(r.split("-")[0]), int(r.split("-")[1])
            range_size = e - s
            if range_size not in obj_dict:
                obj_dict[range_size] = []
            obj_dict[range_size].append({"pdb_id":pdb_id, "s":s, "e":e})
    
    obj_list = []
    for range_size in sorted(obj_dict.keys(), reverse=True):
        obj_list += obj_dict[range_size]


    is_merged = {}
    cls_mat = []
    for i in range(0, len(obj_list)):
        obj_i = obj_list[i]
        obj_i_s = "%s|%s|%s" % (obj_i["pdb_id"], obj_i["s"], obj_i["e"])
        if obj_i_s in is_merged:
            continue
        obj_i["overlap"] = 1.0
        obj_i["overlap_cat"] = 0.75
        cls_list = [obj_i]
        is_merged[obj_i_s] = True
        for j in range(0, len(obj_list)):
            obj_j = obj_list[j]
            obj_j_s = "%s|%s|%s" % (obj_j["pdb_id"], obj_j["s"], obj_j["e"])
            if i == j or obj_j_s in is_merged:
                continue
            min_e, max_s = min(obj_i["e"], obj_j["e"]), max(obj_i["s"], obj_j["s"])
            overlap = round(float(min_e - max_s)/float(obj_i["e"] - obj_i["s"]), 2)
            overlap_cat = 0.00
            overlap_cat = 0.25 if overlap >= 0.25 else overlap_cat
            overlap_cat = 0.50 if overlap >= 0.50 else overlap_cat
            overlap_cat = 0.75 if overlap >= 0.75 else overlap_cat
            if overlap > 0.5 or (obj_j["s"] >= obj_i["s"] and obj_j["e"] <= obj_i["e"]):
                obj_j["overlap"] = overlap
                obj_j["overlap_cat"] = overlap_cat
                cls_list.append(obj_j)
                is_merged[obj_j_s] = True
        cls_mat.append(cls_list)


    return cls_mat

def extract_pdb_map_ds(species):

    fasta_file = "unreviewed/%s_protein_canonicalsequences.fasta" % (species)
    xref_file = "unreviewed/%s_protein_xref_alphafolddb.csv" % (species)

    nt_file = species_obj[species]["nt_file"]
    ebi_dir = "downloads/ebi/current/"
    in_file = "%s%s" % (ebi_dir,nt_file)

    range_dict = {}
    method_dict = {}
    resolution_dict = {}
    with open(in_file, "r") as FR:
        for line in FR:
            if line.find(" <http://purl.uniprot.org/core/chain> ") != -1:
                canon = line.strip().split(" ")[0].split("/")[-1].split("#")[0]
                pdb_id = line.strip().split(" ")[0].split("/")[-1].split("#")[1].split("_")[1].lower()
                chain = line.strip().split(" ")[-2].split("=")[0].replace("\"", "").lower()
                pdb_cmb_id = "%s|%s" % (pdb_id, chain)
                r = line.strip().split(" ")[-2].split("=")[-1].replace("\"","")
                s, e = r.split("-")
                if s.isdigit() and e.isdigit():
                    if int(e) > int(s):
                        if canon not in range_dict:
                            range_dict[canon] = {}
                        if pdb_cmb_id not in range_dict[canon]:
                            range_dict[canon][pdb_cmb_id] = {}
                        range_dict[canon][pdb_cmb_id][r] = True
            elif line.find(" <http://purl.uniprot.org/core/method> ") != -1:
                pdb_id = line.strip().split(" ")[0].split("/")[-1].replace(">", "").lower()
                method = line.split(" ")[-2].split("/")[-1].replace(">", "")
                if pdb_id not in method_dict:
                    method_dict[pdb_id] = {}
                method_dict[pdb_id][method] = True
            elif line.find(" <http://purl.uniprot.org/core/resolution> ") != -1:
                pdb_id = line.strip().split(" ")[0].split("/")[-1].replace(">", "").lower()
                resolution = float(line.strip().split(" ")[2].split("^^")[0].replace("\"", ""))
                if pdb_id not in resolution_dict:
                    resolution_dict[pdb_id] = {}
                resolution_dict[pdb_id][resolution] = True


    method_ordr = ["X-Ray_Crystallography", "NMR_Spectroscopy", "Electron_Microscopy",
        "Infrared_Spectroscopy", "Neutron_Diffraction"
    ]


    row = ["uniprotkb_canonical_ac","sequence_region","pdb_id","pdb_chain","start_pos","end_pos","overlap_ratio"]
    row += ["overlap_category", "experimental_method", "resolution","selection_flag"]
    print "\"%s\"" % ("\",\"".join(row))

    canon_list = list(range_dict.keys())
    #canon_list = ["Q9UKV8-1"]
    for canon in canon_list:
        if canon not in canon_dict:
            continue
        row_dict = {}
        pdb_cmb_id_list = list(range_dict[canon].keys())
        if len(pdb_cmb_id_list) == 1:
            pdb_cmb_id = pdb_cmb_id_list[0]
            pdb_id, chain = pdb_cmb_id.split("|")
            range_list = list(range_dict[canon][pdb_cmb_id].keys())
            method_list = list(method_dict[pdb_id].keys()) if pdb_id in method_dict else []
            min_res = min(list(resolution_dict[pdb_id].keys())) if pdb_id in resolution_dict else 100
            row_list = []
            for method in method_ordr:
                if method not in method_list:
                    continue
                for r in range_list:
                    s, e = r.split("-")
                    region = "region_1"
                    row_list.append([canon,region,pdb_id,chain,str(s),str(e),"1.0","0.75",method, str(min_res)]) 
            for row_idx in range(0, len(row_list)):
                row = row_list[row_idx] + [str(row_idx == 0)]
                s = int(row[4])
                if s not in row_dict:
                    row_dict[s] = []
                row_dict[s].append(row)
        else:
            cls_mat = cluster_pdb_ids(range_dict[canon])
            for i in range(0, len(cls_mat)):
                cls_idx = i + 1
                tmp_dict = {}
                for mem_obj in cls_mat[i]:
                    mem_pdb_id, mem_chain = mem_obj["pdb_id"].split("|")
                    overlap = mem_obj["overlap"]
                    overlap_cat = mem_obj["overlap_cat"]
                    method_list = list(method_dict[mem_pdb_id].keys()) if mem_pdb_id in method_dict else []
                    min_res = min(list(resolution_dict[mem_pdb_id].keys())) if mem_pdb_id in resolution_dict else 100
                    for method in method_list:
                        if method not in tmp_dict:
                            tmp_dict[method] = {}
                        if overlap_cat not in tmp_dict[method]:
                            tmp_dict[method][overlap_cat] = {}
                        if min_res not in tmp_dict[method][overlap_cat]:
                            tmp_dict[method][overlap_cat][min_res] = []
                        tmp_dict[method][overlap_cat][min_res].append(mem_obj)
                row_list = []
                for method in method_ordr:
                    if method not in tmp_dict:
                        continue
                    for overlap_cat in sorted(tmp_dict[method], reverse=True):
                        for min_res in sorted(tmp_dict[method][overlap_cat]):
                            for mem_obj in tmp_dict[method][overlap_cat][min_res]:
                                c, s, e, ovlp = mem_obj["pdb_id"], mem_obj["s"], mem_obj["e"], mem_obj["overlap"]
                                c, chain = c.split("|")
                                region = "region_%s" % (cls_idx)
                                row_list.append([canon,region,c,chain,str(s),str(e),str(ovlp),str(overlap_cat),method, str(min_res)])
                for row_idx in range(0, len(row_list)):
                    row = row_list[row_idx] + [str(row_idx == 0)]
                    s = int(row[4])
                    if s not in row_dict:
                        row_dict[s] = []
                    row_dict[s].append(row)
        for s in sorted(row_dict):
            for row in row_dict[s]:
                print "\"%s\"" % ("\",\"".join(row)) 


    seq_hash = load_fasta_sequences(fasta_file)

    data_frame = {}
    libgly.load_sheet(data_frame, xref_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        xref_id = row[f_list.index("xref_id")]
        ac = canon.split("-")[0]
        if ac != xref_id:
            continue
        pdb_id = "AF-%s-F1" % (canon.split("-")[0])
        pdb_file = "downloads/pdb/current/%s.pdb" % (pdb_id)
        if os.path.isfile(pdb_file) == False:
            continue
        s, e = 1, len(seq_hash[canon])
        ovlp, overlap_cat, method, min_res = 1.0, 1.0, "AlphaFold", -1.0
        region, chain = "global", "global"
        newrow = [canon,region,pdb_id,chain,str(s),str(e),str(ovlp),str(overlap_cat),method, str(min_res), "True"]
        print "\"%s\"" % ("\",\"".join(newrow))

    

    return



def expand_rowcoll(rowcoll, val_list):

    newrowcoll = []
    for row in rowcoll:
        for val in val_list:
            val = val.encode('ascii', 'ignore').decode('ascii') if val != None else ""
            newrowcoll.append(row + [val])
    return newrowcoll



def get_alliance_genome_disease_rowlist(species):
    
    species2xref = {
        "human":"unreviewed/human_protein_xref_hgnc.csv",
        "mouse":"unreviewed/mouse_protein_xref_mgi.csv",
        "rat":"unreviewed/rat_protein_xref_rgd.csv",
        "fruitfly":"unreviewed/fruitfly_protein_xref_flybase.csv",
        "yeast":"unreviewed/yeast_protein_xref_sgd.csv",
        "zebrafish":"unreviewed/zebrafish_protein_xref_zfin.csv"
    }

    row_list = []
    xref2canon = {}
    in_file = species2xref[species]
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        xref2canon[row[f_list.index("xref_id")]] = row[f_list.index("uniprotkb_canonical_ac")]

    seen = {}
    #in_file = "downloads/alliance_genome/current/%s_disease_genome_alliance.tsv" % (species)
    in_file = "downloads/alliance_genome/current/DISEASE-ALLIANCE_COMBINED.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    for row in data_frame["data"]:
        if row[f_list.index("DBobjectType")] != "gene":
            continue
        if len(row) != len(f_list):
            continue
        gene_id = row[f_list.index("DBObjectID")].split(":")[1]
        do_id = row[f_list.index("DOID")].split(":")[1]
        if gene_id not in xref2canon:
            continue
        canon = xref2canon[gene_id]
        combo_id = "%s %s" % (canon, do_id)
        if combo_id not in seen:
            newrow = [canon,do_id]
            row_list.append(newrow)
        seen[combo_id] = True

    return row_list





def get_genomics_england_disease_rowlist(species):

    mimid2doid = {}
    mondoid2mimid = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "protein_disease_idmap.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        mondo_id = row[f_list.index("mondo_id")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "protein_xref_omim":
            mimid2doid[xref_id] = "%s|%s" % (do_id, mondo_id)
            mondoid2mimid[mondo_id] = xref_id


    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_hgnc.csv" % (species)
    if species == "mouse":
        in_file = path_obj["unreviewed"] +  "%s_protein_xref_mgi.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    genename2canon = {}
    genename2hgncid = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_name = row[f_list.index("xref_label")]
        genename2canon[gene_name] = canon
        genename2hgncid[gene_name] = row[f_list.index("xref_id")]


    row_list = [
        ["uniprotkb_canonical_ac","gene_name", "hgnc_id", "mode_of_inheritance",
            "disease_name","pmid","mim_id","doid", "mondo_id"]
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
            d_id = p.split(" ")[-1].split("\t")[-1].replace(")", "")
            if d_id.find(":") != -1:
                prefix, d_id = d_id.split(":")
                if prefix == "MONDO" and d_id in mondoid2mimid:
                    d_id = mondoid2mimid[d_id]

            #print "Robel|%s|%s|%s" % (gene_name, d_id, d_id in mimid2doid)
            if d_id in mimid2doid:
                combo_id = mimid2doid[d_id]
                for pmid in pmid_list:
                    newrow = [canon, gene_name, hgnc_id, inheritance, disease_name, pmid, d_id]
                    newrow += combo_id.split("|")
                    row_list.append(newrow)
    return row_list





def extract_glycosylation_motifs_ds(species):
    
    seq_file = path_obj["unreviewed"] +  "%s_protein_canonicalsequences.fasta" % (species)
    cmd = "readlink -f " + seq_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    prosite_file = path_obj["downloads"] + "/prosite/prosite.dat"
    cmd = "readlink -f " + prosite_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")


    newrow = ["uniprotkb_canonical_ac", "start_pos", "end_pos", "motif"] 
    print "\"%s\"" % ("\",\"".join(newrow))
    
    cmd = "perl /software/ps_scan/ps_scan.pl -d %s -p PS00001 %s " % (prosite_file, seq_file)
    
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
            if str(record.seq.upper()).strip() != "":
                seq_hash[seq_id] = str(record.seq.upper())


    seq_set = {}
    seen_isoform = {}
    for isoform in isoform2canon:
        canon = isoform2canon[isoform]
        if canon not in seq_set:
            seq_str = ">%s\n%s" % (canon, seq_hash[canon])
            seq_set[canon] = [seq_str]
            seen_isoform[canon] = {}
        if isoform != canon and isoform not in seen_isoform[canon]:
            seq_str = ">%s\n%s" % (isoform, seq_hash[isoform])
            seq_set[canon].append(seq_str)
            seen_isoform[canon][isoform] = True


    
    cluster_count = 0
    for canon in seq_set:
        if len(seq_set[canon]) > 1:
            cluster_count += 1

    log_file = "logs/%s_protein_isoform_alignments.log" %(species)
    with open(log_file, "w") as FW:
        FW.write("Started isoform alignments\n")

    
    if os.path.isdir("alignments/isoformset/%s" % (species)) == False:
        cmd = "mkdir alignments/isoformset/%s" % (species)
        x = commands.getoutput(cmd)


    i = 0
    for canon in seq_set:
        if len(seq_set[canon]) > 1:
            seq_file = "alignments/isoformset/%s/isoformset.uniprotkb.%s.fasta" % (species,canon)
            dnd_file = "alignments/isoformset/%s/isoformset.uniprotkb.%s.dnd" % (species,canon)
            aln_file = "alignments/isoformset/%s/isoformset.uniprotkb.%s.aln" % (species,canon)
            seq_buffer_old = ""
            if os.path.isfile(seq_file) == True:
                seq_buffer_old = open(seq_file, "r").read()
            seq_buffer_new = "\n\n".join(seq_set[canon]) + "\n"
            if seq_buffer_new != seq_buffer_old or os.path.isfile(aln_file) == False:
                with open(seq_file, "w") as FW:
                    FW.write("%s" % (seq_buffer_new))
                cmd = "/software/clustalw2/clustalw2 -align -outorder=input -type=protein -infile=%s " % (seq_file)
                x = commands.getoutput(cmd)
                cmd = "rm -f %s" % (dnd_file)
                x = commands.getoutput(cmd)
                i += 1
                with open(log_file, "a") as FL:
                    FL.write("created %s\n" % (aln_file))
            else:
                with open(log_file, "a") as FL: 
                    FL.write("file %s already exists\n" % (aln_file))
        else:
            cmd = "rm -f alignments/isoformset/%s/isoformset.uniprotkb.%s.*" % (species,canon)
            x = commands.getoutput(cmd)


    return







def extract_expression_disease_ds(species):

    work_book = {}
    sheet_name = "masterlist"
    work_book[sheet_name] = {}
    idmap_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")



    in_file = path_obj["downloads"] + "bioxpress/current/bioxpress_disease.csv"
    sheet_name = "diseaseexpression"
    work_book[sheet_name] = {}
    libgly.load_sheet(work_book[sheet_name], in_file, ",")


    row = ["uniprotkb_canonical_ac","significance","direction", "do_id", "do_name", "parent_doid", "parent_doname", "xref_key", "xref_id"]
    print "\"%s\"" % ("\",\"".join(row))

    seen_row = {}
    f_list = work_book[sheet_name]["fields"]
    for row in work_book[sheet_name]["data"]:
        ac = row[f_list.index("uniprotkb_ac")]
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]

        if row[f_list.index("doid")] in ["3963", "0070003", "3119"]:
            continue
        parent_doid, parent_doname = "", ""
        if len(row[f_list.index("parent_name")]) > 0:
            parent_doid, parent_doname = row[f_list.index("parent_doid")], row[f_list.index("parent_name")]
        newrow = [
            canon,
            row[f_list.index("significance")],
            row[f_list.index("direction")],
            row[f_list.index("doid")],
            row[f_list.index("do_name")],
            parent_doid, 
            parent_doname,
            "protein_xref_bioxpress",
            canon.split("-")[0]
        ]
        row_str = json.dumps(newrow)
        if row_str not in seen_row:
            seen_row[row_str] = True
            print "\"%s\"" % ("\",\"".join(newrow))




def extract_expression_normal_ds(species):


    canon2bgee = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_bgee.csv" % (species)
    libgly.load_sheet_as_dict(data_frame, in_file, ",", "uniprotkb_canonical_ac")
    tmp_fl = data_frame["fields"]
    for canon in data_frame["data"]:
        for tmp_row in data_frame["data"][canon][:1]:
            bgee_id = tmp_row[tmp_fl.index("xref_id")]
            canon2bgee[canon] = bgee_id
 

    work_book = {}
    sheet_name = "masterlist"
    work_book[sheet_name] = {}
    idmap_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")
   

    sheet_name = "normalexpression"
    work_book[sheet_name] = {}
    in_file = path_obj["downloads"] + "bioxpress/current/bioxpress_normal.csv"
    libgly.load_sheet(work_book[sheet_name], in_file, ",")        

    selected_fields = [
        "uberon_anatomical_id","uberon_anatomical_name",
        "uberon_developmental_id","uberon_developmental_name", 
        "expression_level_gene_relative","expression_level_anatomical_relative",
        "call_quality","expression_rank_score","expression_score"
    ]
    seen_row = {}
    newrow = ["uniprotkb_canonical_ac"] + selected_fields + ["xref_key", "xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    f_list = work_book[sheet_name]["fields"]
    for row in work_book[sheet_name]["data"]:
        ac = row[f_list.index("uniprotkb_ac")]
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        if canon in canon2bgee:
            bgee_id = canon2bgee[canon]
            newrow = [canon]
            for f in selected_fields:
                newrow.append(row[f_list.index(f)])
            newrow += ["protein_xref_bgee",bgee_id]
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                seen_row[row_str] = True
                print "\"%s\"" % ("\",\"".join(newrow))

    return



def load_ac2canon(species):

    dict_one, dict_two = {}, {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        dict_one[canon] = canon
        isoform_list = [row[f_list.index("reviewed_isoforms")],
                        row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            dict_one[ac] = canon
            dict_one[isoform] = canon
            if isoform == canon:
                dict_two[ac] = canon

    return dict_one, dict_two


def load_canon_dict(species):

    tmp_dict = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        tmp_dict[canon] = True

    return tmp_dict



def load_canon2genename(species):

    canon2genename = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_name = row[f_list.index("gene_name")]
        canon2genename[canon] = gene_name

    return canon2genename





def load_o_glyco_sites(species):

    site_dict = {}
    file_list = glob.glob("unreviewed/%s_proteoform_glycosylation_sites_*.csv" % (species))
    for in_file in file_list:
        if in_file.find("stat") != -1:
            continue
        source = in_file.split("_")[-1].split(".")[0]
        FR = open(in_file, "r")
        idx = 0
        f_list = []
        for line in FR:
            idx += 1
            row = line.strip().split("\",\"")
            row[0] = row[0].replace("\"", "")
            row[-1] = row[-1].replace("\"", "")
            if idx == 1:
                f_list = row
                continue
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            gly_type = row[f_list.index("glycosylation_type")]
            if row[f_list.index("glycosylation_site_uniprotkb")] == "":
                continue
            pos = int(row[f_list.index("glycosylation_site_uniprotkb")].split("|")[0])
            combo_id = "%s %s" % (canon, pos)
            if gly_type.lower() in ["o_linked", "o-linked"]:
                site_dict[combo_id] = True
        FR.close()

    return site_dict

def get_n_site_stat(n_site_info):

    site_stat_dict = {}
    #"glycan":{}, "pmid":{}, "glycancount":{}, "pmidcount":{}}
    for combo_id in n_site_info:
        o = n_site_info[combo_id]
        for o in n_site_info[combo_id]:
            if combo_id not in site_stat_dict:
                site_stat_dict[combo_id]={"glycan":{},"pmid":{},"glycancount":0,"pmidcount":0}
            if o["xrefkey"] == "protein_xref_pubmed":
                site_stat_dict[combo_id]["pmid"][o["xrefid"]] = True
            if o["glytoucanac"] != "":
                site_stat_dict[combo_id]["glycan"][o["glytoucanac"]] = True
        n = len(site_stat_dict[combo_id]["pmid"].keys())
        site_stat_dict[combo_id]["pmidcount"] = n
        
        n = len(site_stat_dict[combo_id]["glycan"].keys())
        site_stat_dict[combo_id]["glycancount"] = n


    return site_stat_dict




def load_canon2genenames(species):

    in_file = path_obj["unreviewed"] + "/%s_protein_genenames_uniprotkb.csv" % (species)
    if os.path.isfile(in_file) == False:
        return {}, {}


    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    canon2genenames = {}
    genename2canon = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        name = row[f_list.index("gene_symbol_recommended")].strip()
        alt = row[f_list.index("gene_symbol_alternative")].strip()
        if name == "":
            continue
        if canon not in canon2genenames:
            canon2genenames[canon] = []
        
        if name not in canon2genenames[canon]:
            canon2genenames[canon].append(name)
        if alt != "" and alt not in canon2genenames[canon]:
            canon2genenames[canon].append(alt)

    in_file = path_obj["unreviewed"] + "/%s_protein_masterlist.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_name = row[f_list.index("gene_name")].strip()
        if gene_name == "":
            continue
        if canon not in canon2genenames:
            canon2genenames[canon] = [gene_name]


    for canon in canon2genenames:
        gene_name = canon2genenames[canon][0]
        genename2canon[gene_name] = canon

    return canon2genenames, genename2canon



def load_canon2recname(species):
    
    canon2recname = {}
    file_list = glob.glob(path_obj["unreviewed"] +  "%s_protein_recnames.csv" % (species))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            canon2recname[canon] = row[f_list.index("recommended_name_full")]

    return canon2recname


def extract_glycogenes_ds(species):

    in_file = "unreviewed/%s_protein_xref_geneid.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    geneid2canon = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_id = row[f_list.index("xref_id")]
        geneid2canon[gene_id] = canon

    in_file = path_obj["downloads"] + "/glycogenes/current/%s_glycogenes.csv" %(species)
    newrow = ["uniprotkb_canonical_ac","gene_id","gene_name","gene_names_alternative",
        "uniprotkb_protein_recommended_name","glycoenzyme_gene_id",
        "gene_group"]
    print "\"%s\"" % ("\",\"".join(newrow))

    seen_row = {}
    FL = open("logs/%s_protein_glycogenes.log" % (species), "w")
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        # Gene Identifier,Gene Symbol,ID,Name,Aliases,Description,Group,
        # Genomic Location,RefSeq,Protein RefSeq,Control?
        glycoenzyme_gene_id = row[f_list.index("Gene Identifier")]
        gene_id_list = row[f_list.index("ID")].strip().split(",")
        gene_group = row[f_list.index("Group")]
        tmp_row_list = []
        if gene_id_list != [""]:
            for gene_id in gene_id_list:
                if gene_id in geneid2canon:
                    canon = geneid2canon[gene_id]
                    tmp_row_list.append([canon,gene_id, canon2genenames[canon][0]])
        if tmp_row_list == []:
            gene_symbol = row[f_list.index("Gene Symbol")]
            if gene_symbol in genename2canon:
                canon = genename2canon[gene_symbol]
                tmp_row_list = [[canon,"", gene_symbol]]
            else:
                FL.write("%s\n" % (",".join(row)))
        
        for tmp_row in tmp_row_list:
            canon = tmp_row[0]
            rec_name = canon2recname[canon] if canon in canon2recname else ""
            alt_gene_names = ""
            if len(canon2genenames[canon]) > 1:
                alt_gene_names = ";".join(canon2genenames[canon][1:])
            newrow = tmp_row + [alt_gene_names,rec_name,glycoenzyme_gene_id,gene_group]
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                seen_row[row_str] = True
                print "\"%s\"" % ("\",\"".join(newrow))
            
    FL.close()

    return



def load_tissueid2name_map():

    in_file = "generated/misc/sample_mapping.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    tmp_dict = {}
    for row in data_frame["data"]:
        if row == []:
            continue
        term_id,term_name = row[f_list.index("term_id")],  row[f_list.index("term_name")]
        sample_type = row[f_list.index("sample_type")]
        tmp_dict[term_id] = term_name

    return tmp_dict



def extract_biomarkers_ds(species):

    tissueid2name = load_tissueid2name_map()
    doid2name = load_doid2name()

    ac2synonyms = {}
    file_list = ["unreviewed/%s_protein_altnames.csv" % (species), "unreviewed/%s_protein_submittednames.csv" % (species)]
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[0].split("-")[0]
            for i in range(1, len(row)):
                name = row[i]
                if ac not in ac2synonyms:
                    ac2synonyms[ac] = []
                if name not in ac2synonyms[ac]:
                    ac2synonyms[ac].append(name)
    

    f_list_one = ["biomarker_id", "biomarker", "assessed_biomarker_entity", "assessed_biomarker_entity_id",
        "assessed_entity_type", "exposure_agent", "exposure_agent_id",
        "best_biomarker_role", "loinc_code", "evidence", "condition"]
    newrow = ["uniprotkb_canonical_ac", "biomarker_canonical_id"]
    newrow += f_list_one
    newrow += ["component_level_tags", "top_level_tags",  "evidence_type", "component_synonyms"]
    newrow += ["uberon_id","anatomical_entity", "do_id","do_name", "do_desc", "do_syn"]
    newrow += ["xref_key","xref_id","src_xref_key", "src_xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))


    log_file = "logs/%s_protein_biomarkers.log" % (species)
    FL = open(log_file, "w")

    #in_file = path_obj["downloads"] + "/biomarkerdb/current/allbiomarkers.tsv"
    in_file = path_obj["downloads"] + "/biomarkerdb/current/oncomx.tsv"



    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        flag_list = []
        biomarker_id = row[f_list.index("biomarker_id")]
        biomarker_canonical_id = row[f_list.index("biomarker_id")].split("-")[0]
        entity_id = row[f_list.index("assessed_biomarker_entity_id")]
        evidence = row[f_list.index("evidence")]
        tag = row[f_list.index("tag")]
        evidence_source = row[f_list.index("evidence_source")]
        do_id = row[f_list.index("condition_id")].lower().replace("doid:", "")
        uberon_id = row[f_list.index("specimen_id")].lower().replace("uberon:", "")
        if entity_id.find("UPKB:") == -1:
            flag_list.append("no-UPKB-component")
        ac = entity_id.split(":")[1]
        if ac not in ac2canon:
            flag_list.append("no-ac2canon-map")
        k =  evidence_source.split(":")[0].lower()
        if k not in ["pubmed", "doi"]:
            flag_list.append("no-pubmed/doi-xref")

        if flag_list != []:
            FL.write("\"%s\"\n" % ("\",\"".join(row + [";".join(flag_list)])) )
            continue


        synonyms = ac2synonyms[ac] if ac in ac2synonyms else []
        component_synonyms = "|".join(synonyms)
        canon = ac2canon[ac]
        xref_key = "protein_xref_" + k
        xref_id = evidence_source.split(":")[1]
        tmp_row = [canon, biomarker_canonical_id]
        tag_list_one, tag_list_two = [], []
        for f in f_list_one:
            val = "" if f == "loinc_code" and val == "NA" else row[f_list.index(f)]
            tmp_row.append(val)
        for t in tag.split(";"):
            if t in ["condition","best_biomarker_role"]:
                tag_list_two.append(t)
            else:
                tag_list_one.append(t)
        n_one = len(tag_list_one)
        n_two = len(tag_list_two)
        evidence_type = ""
        if n_one > 0 and n_two > 0:
            evidence_type = "both_levels"
        elif n_one > 0:
            evidence_type = "component_level"
        elif n_two > 0:
            evidence_type = "top_level"
        tmp_row.append("; ".join(tag_list_one))
        tmp_row.append("; ".join(tag_list_two))
        tmp_row.append(evidence_type)
        tmp_row.append(component_synonyms)

        uberon_name, do_name  = "", ""
        if do_id in doid2name["do"]:
            do_name = doid2name["do"][do_id]["name"]
       
        do_desc = doid2name["do"][do_id]["description"] if do_id in doid2name["do"] else ""
        do_syn_list = doid2name["do"][do_id]["synonyms"] if do_id in doid2name["do"] else [] 
        do_syn = "|".join(do_syn_list)
        u_id = "UBERON:" + uberon_id 
        if u_id in tissueid2name:
            uberon_name = tissueid2name[u_id]

        tmp_row += [uberon_id,uberon_name,do_id,do_name, do_desc, do_syn]
        #src_xref_key = "protein_xref_glygen_ds"
        #ds = "%s_protein_biomarkers" % (species)
        #src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        src_xref_key = "protein_xref_biomarkerkb"
        src_xref_id = biomarker_id 
        newrow = tmp_row + [xref_key, xref_id, src_xref_key, src_xref_id]
        if xref_id != "https":
            print "\"%s\"" % ("\",\"".join(newrow))

    FL.close()

    return


def extract_mutation_germline_ds(species):

    fasta_file = "unreviewed/%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(fasta_file)
   
    n_site_info = get_glycosylation_sites(species)
    is_o_site = load_o_glyco_sites(species)
    site_stat_dict = get_n_site_stat(n_site_info)


    is_glycoprotein = {}
    for combo_id in n_site_info.keys() + is_o_site.keys():
        canon = combo_id.split(" ")[0]
        is_glycoprotein[canon] = True

    data_frame = {}
    in_file = "unreviewed/protein_disease_idmap.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    mimid2doid = {}
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "protein_xref_omim":
            if xref_id not in mimid2doid:
                mimid2doid[xref_id] = []
            mimid2doid[xref_id].append(do_id)
   
    tax_name = species_obj[species]["long_name"].lower().replace(" ", "-")
    in_file = path_obj["downloads"] + "/ebi/current/dbSNP-%s.tsv" % (tax_name)

    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    seen_row = {}
    import sys
    reload(sys)
    sys.setdefaultencoding('utf-8')
    import io
    with io.open(in_file, "r", encoding="utf-8",errors="ignore") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            row = line.replace("\"", "\'").split("\t")
            if lcount == 1:
                f_list = row
                header_row = ["uniprotkb_canonical_ac", "aa_pos"]
                for f in f_list:
                    if f.strip() not in ["uniprotkb_accession","data_source",""]:
                        f = "ref_nt" if f.strip() == "ref_allele" else f
                        f = "alt_nt" if f.strip() == "alt_allele" else f
                        f = "chr_id" if f.strip() == "chromosome_id" else f
                        f = "chr_pos" if f.strip() == "position" else f
                        f = "minor_allelic_frequency" if f.strip() == "frequency" else f
                        header_row.append(f.strip())
                header_row += ["xref_key", "xref_id", "mim_id", "do_id", "filter_flags"]
                header_row += ["glyco_annotation"]
                print "\"%s\"" % ("\",\"".join(header_row))
            else:
                if len(row) != len(f_list):
                    continue
                if row[f_list.index("somatic_status")] == "1":
                    continue
                ac = row[f_list.index("uniprotkb_accession")]
                if ac.find("-") != -1 and ac not in canon_dict:
                    continue
                if ac not in ac2canon:
                    continue
                canon = ac2canon[ac]
                dbsnp_id = row[f_list.index("dbsnp_id")]
                mim_id_list = row[f_list.index("disease_xrefs")].strip().split(",")
                begin_aa_pos = int(row[f_list.index("begin_aa_pos")])
                end_aa_pos = int(row[f_list.index("end_aa_pos")])
                ref_aa = row[f_list.index("ref_aa")]
                alt_aa = row[f_list.index("alt_aa")]
                aa_pos = begin_aa_pos
                
                #consider only point mutations
                if begin_aa_pos != end_aa_pos:
                    continue
                if begin_aa_pos > len(seq_hash[canon]):
                    continue
                #check amino acid
                if ref_aa != seq_hash[canon][begin_aa_pos-1]:
                    continue
 
                
                evdn_list = row[f_list.index("evidence_ECO:0000313")].split(",")
                data_source_list = [] 
                row_one =  [canon, str(aa_pos)]
                for f in f_list:
                    if f.strip() == "data_source":
                        data_source_list = row[f_list.index(f.strip())].strip().split(",")
                    if f.strip() in ["uniprotkb_accession", "data_source",""]:
                        continue
                    f_val = row[f_list.index(f.strip())].strip()
                    if f == "frequency" and f_val.find("E-") != -1:
                        f_val = str(round(float(f_val), 6))
                    row_one.append(f_val)

                map_list = []
                if mim_id_list != [""]:
                    for idx in mim_id_list:
                        mim_id = idx.split(":")[1].strip()
                        if mim_id in mimid2doid:
                            for do_id in mimid2doid[mim_id]:
                                map_list.append([mim_id, do_id])
                

                e_list = get_mutation_effect_list(seq_hash[canon], aa_pos,ref_aa, alt_aa)
                
                filter_flag_list = []
                for o in e_list:
                    o_combo = "%s %s" % (canon, aa_pos)
                    ref_m,alt_m,eff = o["refmotif"], o["altmotif"], o["effect"]
                    if eff == "o-glyco-site-loss" and o_combo not in is_o_site:
                        continue
                    n_combo = "%s %s" % (canon, o["motifstart"])
                    if eff == "n-glyco-sequon-loss" and n_combo not in n_site_info:
                        continue
                    if eff == "n-glyco-sequon-gain" and canon not in is_glycoprotein:
                        continue
                    filter_flag_list.append("%s (%s->%s)" % (eff, ref_m,alt_m))
                
                if filter_flag_list == []:
                    continue

                filter_flags = ""
                filter_flags = "Germline mutation passed %s out of %s filters: " % (1, 1)
                filter_flags += "; ".join(filter_flag_list) + "."
                
                site_id = "%s %s" % (canon, aa_pos)
                glycan_count, pmid_count = 0, 0
                if site_id in site_stat_dict:
                    glycan_count = site_stat_dict[site_id]["glycancount"]
                    pmid_count = site_stat_dict[site_id]["pmidcount"]
                glyco_annotation = "glycancount:%s, pmidcount:%s" % (glycan_count,pmid_count)
                if map_list != []:
                    for map_row in map_list:
                        xref_key, xref_id = "protein_xref_dbsnp", dbsnp_id
                        mim_id, do_id =  map_row[0], map_row[1]
                        row_two = [xref_key, xref_id, mim_id,do_id]
                        row_two += [filter_flags,glyco_annotation]
                        row_str = json.dumps(row_one + row_two)
                        if row_str not in seen_row:
                            seen_row[row_str] = True
                            print "\"%s\"" % ("\",\"".join(row_one + row_two))
                        for dsrc in data_source_list:
                            xref_key = "protein_xref_" + dsrc.lower()
                            xref_id = dbsnp_id
                            row_two = [xref_key, xref_id, mim_id,do_id]
                            row_two += [filter_flags,glyco_annotation]
                            row_str = json.dumps(row_one + row_two)
                            if row_str not in seen_row:
                                seen_row[row_str] = True
                                print "\"%s\"" % ("\",\"".join(row_one + row_two))
                        for evdn in evdn_list:
                            if evdn.find("pubmed") == -1:
                                continue
                            xref_key = "protein_xref_pubmed"
                            xref_id = evdn.split(":")[1]
                            row_two = [xref_key, xref_id, mim_id,do_id]
                            row_two += [filter_flags,glyco_annotation]
                            row_str = json.dumps(row_one + row_two)
                            if row_str not in seen_row:
                                seen_row[row_str] = True
                                print "\"%s\"" % ("\",\"".join(row_one + row_two))
                elif filter_flags != "":
                    xref_key, xref_id = "protein_xref_dbsnp", dbsnp_id
                    mim_id, do_id = "", ""
                    row_two = [xref_key, xref_id, mim_id,do_id]
                    row_two += [filter_flags,glyco_annotation]
                    row_str = json.dumps(row_one + row_two)
                    if row_str not in seen_row:
                        seen_row[row_str] = True
                        print "\"%s\"" % ("\",\"".join(row_one + row_two))
                    for dsrc in data_source_list:
                        xref_key = "protein_xref_" + dsrc.lower()
                        xref_id = dbsnp_id
                        row_two = [xref_key, xref_id, mim_id,do_id] 
                        row_two += [filter_flags,glyco_annotation]
                        row_str = json.dumps(row_one + row_two)
                        if row_str not in seen_row:
                            seen_row[row_str] = True
                            print "\"%s\"" % ("\",\"".join(row_one + row_two))
                    for evdn in evdn_list:
                        if evdn.find("pubmed") == -1:
                            continue
                        xref_key = "protein_xref_pubmed"
                        xref_id = evdn.split(":")[1]
                        row_two = [xref_key, xref_id, mim_id,do_id]
                        row_two += [filter_flags,glyco_annotation]
                        row_str = json.dumps(row_one + row_two)
                        if row_str not in seen_row:
                            seen_row[row_str] = True
                            print "\"%s\"" % ("\",\"".join(row_one + row_two))

    return

def extract_mutation_literature_ds(species):

       
    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_canonicalsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())


    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_dict = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        one = row[f_list.index("one")]
        three = row[f_list.index("three")]
        aa_dict[three.upper()] = one.upper()
        aa_dict[one.upper()] = one.upper()

    data_frame = {}
    #in_file = path_obj["downloads"] +  "biomuta/current/human_cancer_mutation_literature.csv"
    in_file = "compiled/human_cancer_mutation_literature.csv"
    
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    newrow = ["uniprotkb_canonical_ac"]
    for f in f_list:
        if f not in ["pmid","uniprotkb_ac", "gene_symbol", "mutation_mention"]:
            newrow.append(f)
    newrow += ["start_pos", "start_pos", "ref_aa", "alt_aa", "xref_key","xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))

    seen_row = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("uniprotkb_ac")]
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        newrow = [canon]
        for f in f_list:
            if f not in ["pmid", "uniprotkb_ac", "gene_symbol", "mutation_mention"]:
                newrow.append(row[f_list.index(f)])
        xref_key, xref_id = "protein_xref_pubmed", row[f_list.index("pmid")]

        mutation_mention = row[f_list.index("mutation_mention")]
        parts = re.split(r'[0-9]+', mutation_mention)
        pos = mutation_mention.replace(parts[0], "")
        pos = pos.replace(parts[1], "")
        ref_aa = parts[0].upper()
        alt_aa = parts[1].upper()
        if ref_aa not in aa_dict or alt_aa not in aa_dict:
            continue
        start_pos = int(pos)
        end_pos = start_pos
        ref_aa = aa_dict[ref_aa]
        alt_aa = aa_dict[alt_aa]
        if start_pos > len(seq_hash[canon]):
            continue
        if seq_hash[canon][start_pos-1] != ref_aa:
            continue
        newrow += [str(start_pos), str(end_pos), ref_aa, alt_aa, xref_key, xref_id]
        row_str = json.dumps(newrow)
        if row_str not in seen_row:
            seen_row[row_str] = True
            print "\"%s\"" % ("\",\"".join(newrow))
    return


def load_fasta_sequences(fasta_file):
    seq_hash = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())
    return seq_hash




def extract_mutation_somatic_ds(species):

    fasta_file = "unreviewed/%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(fasta_file)

    n_site_info = get_glycosylation_sites(species) 
    is_o_site = load_o_glyco_sites(species)
   
    site_stat_dict = get_n_site_stat(n_site_info)

    is_glycoprotein = {}
    for combo_id in n_site_info.keys() + is_o_site.keys():
        canon = combo_id.split(" ")[0]
        is_glycoprotein[canon] = True
    



    seen_row = {}
    data_frame = {}
    in_file = "unreviewed/protein_disease_idmap.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    mimid2doid = {}
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "protein_xref_omim":
            if xref_id not in mimid2doid:
                mimid2doid[xref_id] = []
            mimid2doid[xref_id].append(do_id)
    

    clinvar_ann = {}
    tax_name = species_obj[species]["long_name"].lower().replace(" ", "-")
    in_file = path_obj["downloads"] + "/ebi/current/dbSNP-%s.tsv" % (tax_name)
    

 
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    import sys
    reload(sys)
    sys.setdefaultencoding('utf-8')
    import io
    with io.open(in_file, "r", encoding="utf-8",errors="ignore") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            row = line.replace("\"", "\'").split("\t")
            if lcount == 1:
                f_list = row
                header_row = ["uniprotkb_canonical_ac", "aa_pos"]
                for f in f_list:
                    if f.strip() not in ["uniprotkb_accession", "data_source",""]:
                        f = "ref_nt" if f.strip() == "ref_allele" else f
                        f = "alt_nt" if f.strip() == "alt_allele" else f
                        f = "chr_id" if f.strip() == "chromosome_id" else f
                        f = "chr_pos" if f.strip() == "position" else f
                        f = "minor_allelic_frequency" if f.strip() == "frequency" else f 
                        header_row.append(f.strip())
                header_row += ["xref_key", "xref_id", "mim_id", "do_id", "filter_flags"]
                header_row += ["glyco_annotation"]
                print "\"%s\"" % ("\",\"".join(header_row))
            else:
                if len(row) != len(f_list):
                    continue
                if row[f_list.index("somatic_status")] == "0":
                    continue
                ac = row[f_list.index("uniprotkb_accession")]
                if ac.find("-") != -1 and ac not in canon_dict:
                    continue
                if ac not in ac2canon:
                    continue
                canon = ac2canon[ac]
                mim_id_list = row[f_list.index("disease_xrefs")].strip().split(",")
                begin_aa_pos = int(row[f_list.index("begin_aa_pos")])
                end_aa_pos = int(row[f_list.index("end_aa_pos")])
                ref_aa = row[f_list.index("ref_aa")]
                alt_aa = row[f_list.index("alt_aa")]
                aa_pos = begin_aa_pos
                #consider only point mutations
                if begin_aa_pos != end_aa_pos:
                    continue
                if begin_aa_pos > len(seq_hash[canon]):
                    continue
                #check amino acid
                if ref_aa != seq_hash[canon][begin_aa_pos-1]:
                    continue
                evdn_list = row[f_list.index("evidence_ECO:0000313")].split(",")
                data_source_list = []
                row_one =  [canon, str(aa_pos)]
                for f in f_list:
                    if f.lower() == "data_source":
                        data_source_list = row[f_list.index(f.strip())].strip().split(",")
                    if f.strip() in ["uniprotkb_accession", "", "data_source"]:
                        continue
                    f_val = row[f_list.index(f.strip())].strip()
                    if f == "frequency" and f_val.find("E-") != -1:
                        f_val = str(round(float(f_val), 6))
                    row_one.append(f_val)
                dbsnp_id = row[f_list.index("dbsnp_id")]
                e_list = get_mutation_effect_list(seq_hash[canon], aa_pos,ref_aa, alt_aa)
                filter_flag_list = []
                for o in e_list:
                    o_combo = "%s %s" % (canon, aa_pos)
                    ref_m,alt_m,eff = o["refmotif"], o["altmotif"], o["effect"]
                    if eff == "o-glyco-site-loss" and o_combo not in is_o_site:
                        continue
                    n_combo = "%s %s" % (canon, o["motifstart"])
                    if eff == "n-glyco-sequon-loss" and n_combo not in n_site_info:
                        continue
                    if eff == "n-glyco-sequon-gain" and canon not in is_glycoprotein:
                        continue
                    filter_flag_list.append("%s (%s->%s)" % (eff,ref_m,alt_m))
                if filter_flag_list == []:
                    continue
                filter_flags = ""
                filter_flags = "Somatic mutation passed %s out of %s filters: " % (1, 1)
                filter_flags += "; ".join(filter_flag_list) + "."
               
                site_id = "%s %s" % (canon, aa_pos)
                glycan_count, pmid_count = 0, 0
                if site_id in site_stat_dict:
                    glycan_count = site_stat_dict[site_id]["glycancount"]
                    pmid_count = site_stat_dict[site_id]["pmidcount"]
                glyco_annotation = "glycancount:%s, pmidcount:%s" % (glycan_count,pmid_count)
                

                if mim_id_list != [""]:
                    for idx in mim_id_list:
                        mim_id = idx.split(":")[1].strip()
                        if mim_id in mimid2doid:
                            for do_id in mimid2doid[mim_id]:
                                xref_key, xref_id = "protein_xref_dbsnp", dbsnp_id
                                row_two = [xref_key, xref_id, mim_id,do_id]
                                row_two += [filter_flags,glyco_annotation]
                                row_str = json.dumps(row_one + row_two)
                                if row_str not in seen_row:
                                    seen_row[row_str] = True
                                    print "\"%s\"" % ("\",\"".join(row_one + row_two))
                                for dsrc in data_source_list:
                                    xref_key = "protein_xref_" + dsrc.lower()
                                    xref_id = dbsnp_id
                                    row_two = [xref_key, xref_id, mim_id,do_id]
                                    row_two += [filter_flags,glyco_annotation]
                                    row_str = json.dumps(row_one + row_two)
                                    if row_str not in seen_row:
                                        seen_row[row_str] = True
                                        print "\"%s\"" % ("\",\"".join(row_one + row_two))
                                for evdn in evdn_list:
                                    if evdn.find("pubmed") == -1:
                                        continue
                                    xref_key = "protein_xref_pubmed"
                                    xref_id = evdn.split(":")[1]
                                    row_two = [xref_key, xref_id, mim_id,do_id]
                                    row_two += [filter_flags,glyco_annotation]
                                    row_str = json.dumps(row_one + row_two)
                                    if row_str not in seen_row:
                                        seen_row[row_str] = True
                                        print "\"%s\"" % ("\",\"".join(row_one + row_two))
                elif filter_flags != "":
                    xref_key, xref_id = "protein_xref_dbsnp", dbsnp_id
                    mim_id, do_id = "", ""
                    row_two = [xref_key, xref_id, mim_id,do_id]
                    row_two += [filter_flags,glyco_annotation]
                    row_str = json.dumps(row_one + row_two)
                    if row_str not in seen_row:
                        seen_row[row_str] = True 
                        print "\"%s\"" % ("\",\"".join(row_one + row_two))
                    for dsrc in data_source_list:
                        xref_key = "protein_xref_" + dsrc.lower()
                        xref_id = dbsnp_id
                        row_two = [xref_key, xref_id, mim_id,do_id]
                        row_two += [filter_flags,glyco_annotation]
                        row_str = json.dumps(row_one + row_two)
                        if row_str not in seen_row:
                            seen_row[row_str] = True
                            print "\"%s\"" % ("\",\"".join(row_one + row_two))
                    for evdn in evdn_list:
                        if evdn.find("pubmed") == -1:
                            continue
                        xref_key = "protein_xref_pubmed"
                        xref_id = evdn.split(":")[1]
                        row_two = [xref_key, xref_id, mim_id,do_id]
                        row_two += [filter_flags,glyco_annotation]
                        row_str = json.dumps(row_one + row_two)
                        if row_str not in seen_row:
                            seen_row[row_str] = True
                            print "\"%s\"" % ("\",\"".join(row_one + row_two))


    return


def extract_mutation_cancer_ds(species):

    
    fasta_file = "unreviewed/%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(fasta_file)
   
    dbsnp_file = path_obj["downloads"] + "/ebi/current/toy.tsv"
    biomuta_file = path_obj["downloads"] + "/biomuta/current/toy.csv"
    log_file = "logs/%s_protein_mutation_cancer.log" % (species)
    clinical_file = path_obj["downloads"] + "/biomuta/current/clinical_information.csv"
    n_site_info, is_o_site, site_stat_dict, is_glycoprotein = {}, {}, {}, {}
    if DEBUG == False:
        tax_name = species_obj[species]["long_name"].lower().replace(" ", "-")
        dbsnp_file = path_obj["downloads"] + "/ebi/current/dbSNP-%s.tsv" % (tax_name)
        biomuta_file = path_obj["downloads"] + "/biomuta/current/biomuta.csv"
        n_site_info = get_glycosylation_sites(species)

        is_o_site = load_o_glyco_sites(species)
        site_stat_dict = get_n_site_stat(n_site_info)
        is_glycoprotein = {}
        for combo_id in n_site_info.keys() + is_o_site.keys():
            canon = combo_id.split(" ")[0]
            is_glycoprotein[canon] = True

    seen_canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        seen_canon[canon] = True

    in_file = "unreviewed/protein_disease_idmap.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    
    mimid2doid = {}
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "protein_xref_omim":
            if xref_id not in mimid2doid:
                mimid2doid[xref_id] = []
            mimid2doid[xref_id].append(do_id)

    lit_ann = {}
    lit_list = {}
    in_file = "unreviewed/human_protein_mutation_literature.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        ref_aa = row[f_list.index("ref_aa")]
        alt_aa = row[f_list.index("alt_aa")]
        aa_pos = row[f_list.index("start_pos")]
        do_id = row[f_list.index("doid")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "protein_xref_pubmed":
            pmid = xref_id
            combo_id = "%s %s %s %s" % (canon, aa_pos, ref_aa, alt_aa)
            if combo_id not in lit_ann:
                lit_ann[combo_id] = []
                lit_list[combo_id] = []
            lit_ann[combo_id].append(do_id)
            if pmid.strip() != "" and pmid not in lit_list[combo_id]:
                lit_list[combo_id].append(pmid)


    clinical_info = {}
    cmd = "readlink -f " + clinical_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")
    #These info come from TCGA. The file that needs to be parsed is glygen/downloads/biomuta/v-5.0/tcga/clinical_information.csv. Target the fields submitter_id for patient and project_id for cancer. The submitter id is present in compiled/biomuta_v5.csv as sample_name to map the info.
    libgly.load_sheet(data_frame, clinical_file, ",")
    f_list = data_frame["fields"]
    studyid2tested = {}
    patientid2studyid = {}
    for row in data_frame["data"]:
        patient_id = row[f_list.index("submitter_id")]
        study_id = row[f_list.index("proj__project_id")]
        if study_id not in studyid2tested:
            studyid2tested[study_id] = 0
        studyid2tested[study_id] += 1
        patientid2studyid[patient_id] = study_id
    

    
    
    clinvar_ann = {}
    dbsnpid_dict = {}
    cmd = "readlink -f " + dbsnp_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")


    seen_row = {}
    import sys
    reload(sys)
    sys.setdefaultencoding('utf-8')
    import io
    with io.open(dbsnp_file, "r", encoding="utf-8",errors="ignore") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            row = line.split("\t")
            if lcount == 1:
                f_list = row
            else:
                if row[f_list.index("somatic_status")] == "0":
                    continue
                ac = row[f_list.index("uniprotkb_accession")].split("-")[0]
                if ac not in ac2canon:
                    continue
                canon = ac2canon[ac]
                aa_pos = int(row[f_list.index("begin_aa_pos")])
                ref_aa = row[f_list.index("ref_aa")]
                alt_aa = row[f_list.index("alt_aa")]
                mim_id_list = row[f_list.index("disease_xrefs")].strip().split(",")
                dbsnp_id = row[f_list.index("dbsnp_id")]
                combo_id = "%s %s %s %s" % (canon, aa_pos, ref_aa, alt_aa)
                
                #we get dbSNP rs ID from Biomuta file now
                #if dbsnp_id.strip() != "":
                #    if combo_id not in dbsnpid_dict:
                #        dbsnpid_dict[combo_id] = {}
                #    dbsnpid_dict[combo_id][dbsnp_id] = True
                if mim_id_list != [""]:
                    for idx in mim_id_list:
                        mim_id = idx.split(":")[1].strip()
                        if mim_id in mimid2doid:
                            if combo_id not in clinvar_ann:
                                clinvar_ann[combo_id] = []
                            clinvar_ann[combo_id] += mimid2doid[mim_id]



    FL = open(log_file, "w")

    row = ["uniprotkb_canonical_ac","aa_pos","ref_aa","alt_aa","chr_id","chr_start_pos",
            "chr_end_pos",
            "ref_nt","alt_nt", "patients_positive", "patients_tested", "mut_freq",
            "data_source","do_id","do_name","xref_key", "xref_id", "filter_flags", 
            "glyco_annotation", "mim_id", "somatic_status","minor_allelic_frequency"
    ]
    print "\"%s\"" % ("\",\"".join(row))

    anchor_fields = ["uniprotkb_canonical_ac","aa_pos","ref_aa","alt_aa","do_id","source"]
    cmd = "readlink -f " + biomuta_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    f_list = []
    cancer_count_dict = {}
    source_dict = {}
    data_grid = {}
    mut2patient = {}
    
    with io.open(biomuta_file, "r", encoding="utf-8",errors="ignore") as FR:
        row_count = 0
        for line in FR:
            row = line.strip().split("\",\"")
            row[0] = row[0].replace("\"", "")
            row[-1] = row[-1].replace("\"", "")
    #with open(biomuta_file, "r") as FR:
    #    csv_grid = csv.reader(FR, delimiter=",", quotechar='\"')
    #    row_count = 0
    #    for row in csv_grid:
            row_count += 1
            if row_count == 1:
                f_list = row
                newrow = row + ["flag_list"]
                FL.write("\"%s\"\n" % ("\",\"".join(newrow)))
                continue
            flag_list = []
            sample_name = row[f_list.index("sample_name")]
            chr_id = row[f_list.index("chr_id")]
            start_pos = row[f_list.index("start_pos")]
            end_pos = row[f_list.index("end_pos")]
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            aa_pos = row[f_list.index("aa_pos")]
            ref_nt = row[f_list.index("ref_nt")]
            alt_nt = row[f_list.index("alt_nt")]
            ref_aa = row[f_list.index("ref_aa")]
            alt_aa = row[f_list.index("alt_aa")]
            do_name = row[f_list.index("do_name")]
            source = row[f_list.index("source")]
            if source == "cosmic":  
                continue

            if row[f_list.index("dbsnp_id")].strip() != "":
                dbsnp_id_list = row[f_list.index("dbsnp_id")].split("|")
                for dbsnp_id in dbsnp_id_list:
                    if dbsnp_id[:2] == "rs":
                        combo_id = "%s %s %s %s" % (canon, aa_pos, ref_aa, alt_aa)
                        if combo_id not in dbsnpid_dict:
                            dbsnpid_dict[combo_id] = {}
                        dbsnpid_dict[combo_id][dbsnp_id] = True


            
            
            #do_id = row[f_list.index("do_id")]
            do_id = ""
            if do_name.find("DOID:") != -1:
                do_id = do_name.split("/")[0].split(":")[1].strip()

            if do_id == "":
                flag_list.append("no-disease-info")
            if start_pos.isdigit() == False:
                flag_list.append("start_pos-bad-value")
            if end_pos.isdigit() == False:
                flag_list.append("end_pos-bad-value")

            if canon not in seen_canon:
                flag_list.append("canon-not-found")
            elif canon not in seq_hash:
                flag_list.append("no-canon-seq")
            else:
                if aa_pos.isdigit() == False:
                    flag_list.append("aa_pos-bad-value")
                else:
                    if int(aa_pos) < 1 or int(aa_pos) > len(seq_hash[canon]):
                        flag_list.append("aa_pos > canon_seq_len")
                    elif ref_aa != seq_hash[canon][int(aa_pos)-1]:
                        flag_list.append("ref_aa-mismatch")

            if ref_nt == alt_nt:
                flag_list.append("no-nt-change")
            if ref_aa == alt_aa:
                flag_list.append("no-aa-change")

            if flag_list != []:
                newrow = row + [";".join(flag_list)]
                FL.write("\"%s\"\n" % ("\",\"".join(newrow)))
                continue

            combo_id = "%s %s %s %s" % (canon, aa_pos, ref_aa, alt_aa)
            if combo_id not in cancer_count_dict:
                cancer_count_dict[combo_id] = []
            if do_id not in cancer_count_dict[combo_id]:
                cancer_count_dict[combo_id].append(do_id)
            
            src_combo_id = "%s %s %s %s %s" % (canon,aa_pos,ref_aa,alt_aa,do_id)
            if src_combo_id not in source_dict:
                source_dict[src_combo_id] = []
            if source not in source_dict[src_combo_id]:
                source_dict[src_combo_id].append(source)
            
            new_combo_id = "%s %s %s %s %s %s" % (canon,aa_pos,ref_aa,alt_aa,do_id,source)
            
            if new_combo_id not in data_grid:
                data_grid[new_combo_id] = {}
            for f in f_list:
                if f not in anchor_fields:
                    data_grid[new_combo_id][f] = row[f_list.index(f)]
            if new_combo_id not in mut2patient:
                mut2patient[new_combo_id] = {}
            mut2patient[new_combo_id][sample_name] = True


    FL.close()


    freq_dist = {}
    #cutoff_one = 0.01
    #cutoff_two = 0
    cutoff_one = 1.0
    cutoff_two = 10
    for new_combo_id in data_grid:
        obj = data_grid[new_combo_id]
        canon,aa_pos, ref_aa, alt_aa,do_id,source = new_combo_id.split(" ")
        do_name = obj["do_name"]
        patient_id = data_grid[new_combo_id]["sample_name"]
        patients_tested, patients_positive, mut_freq = 0, 0, 0.0
        if patient_id in patientid2studyid:
            study_id = patientid2studyid[patient_id]
            patient_list = list(mut2patient[new_combo_id].keys())
            patients_tested = studyid2tested[study_id] 
            patients_positive = len(patient_list)
            mut_freq = round(100.0*float(patients_positive)/float(patients_tested), 1)
        
        if ref_aa == alt_aa:
            continue
        combo_id = "%s %s %s %s" % (canon, aa_pos, ref_aa, alt_aa)
        src_combo_id = "%s %s %s %s %s" % (canon, aa_pos, ref_aa, alt_aa,do_id)

        if src_combo_id in source_dict:
            if source == "icgc" and "tcga" in source_dict[src_combo_id]:
                continue

        cancer_count = len(cancer_count_dict[combo_id])
        filter_flag_list = []
        if mut_freq >= cutoff_one:
            flag_label = "patient freq. (%s%s)" % (mut_freq, "%")
            filter_flag_list.append(flag_label)
        if patients_positive >= cutoff_two:
            flag_label = "patient count (%s/%s)"%(patients_positive,patients_tested)
            filter_flag_list.append(flag_label)
        if cancer_count >= 3:
            flag_label = "num. of cancers (%s)" % (cancer_count)
            filter_flag_list.append(flag_label)
        if combo_id in clinvar_ann:
            if do_id in clinvar_ann[combo_id]:
                flag_label = "DB annot. (%s)" % (1)
                filter_flag_list.append(flag_label)
        if combo_id in lit_ann:
            if do_id in lit_ann[combo_id]:
                flag_label = "Lit. mined (%s)" % (len(lit_list[combo_id]))
                filter_flag_list.append(flag_label)


        e_list = get_mutation_effect_list(seq_hash[canon], int(aa_pos),ref_aa, alt_aa)
        extra_filter_flag_list = []
        for o in e_list:
            o_combo = "%s %s" % (canon, aa_pos)
            ref_m,alt_m,eff = o["refmotif"],o["altmotif"],o["effect"]
            if eff == "o-glyco-site-loss" and o_combo not in is_o_site:
                continue
            n_combo = "%s %s" % (canon, o["motifstart"])
            if eff == "n-glyco-sequon-loss" and n_combo not in n_site_info:
                continue
            if eff == "n-glyco-sequon-gain" and canon not in is_glycoprotein:
                continue
            extra_filter_flag_list.append("%s (%s->%s)" % (eff,ref_m,alt_m))

        if filter_flag_list + extra_filter_flag_list == []:
            continue
      
        n_filters_passed = len(filter_flag_list)
        if extra_filter_flag_list != []:
            n_filters_passed += 1

        filter_flags = "Somatic mutation passed %s out of %s filters: " % (n_filters_passed, 6)
        filter_flags += "; ".join(filter_flag_list + extra_filter_flag_list) + "."
        if combo_id in dbsnpid_dict:
            rs_list = list(dbsnpid_dict[combo_id].keys())
            if rs_list != []:
                filter_flags += " In dbSNP: %s." % ("; ".join(rs_list))

        if canon not in seen_canon:
            continue
        if do_id in ["3963", "0070003", "3119"]:
            continue

        site_id = "%s %s" % (canon, aa_pos)
        glycan_count, pmid_count = 0, 0
        if site_id in site_stat_dict:
            glycan_count = site_stat_dict[site_id]["glycancount"]
            pmid_count = site_stat_dict[site_id]["pmidcount"]
        glyco_annotation = "glycancount:%s, pmidcount:%s" % (glycan_count,pmid_count)
        newrow = [
            canon,
            str(aa_pos),
            ref_aa,
            alt_aa,
            obj["chr_id"],
            obj["start_pos"],
            obj["end_pos"],
            obj["ref_nt"],
            obj["alt_nt"],
            str(patients_positive),
            str(patients_tested),
            str(mut_freq),
            source,
            do_id,
            do_name,
            "protein_xref_biomuta",
            canon.split("-")[0],
            filter_flags,
            glyco_annotation,
            "",
            "1",
            ""
        ]
        row_str = json.dumps(newrow)
        if row_str not in seen_row:
            seen_row[row_str] = True
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


def extract_matrixdb_ds(species):


    glycan_dict = load_glycan_masterdict() 

    in_file = "generated/misc/matrix_db_label.csv"
    matrixdbid2label = {}
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        matrix_db_id = row[f_list.index("matrix_db_id")]
        matrixdbid2label[matrix_db_id] = row[f_list.index("matrix_db_label")]


    matrixdbid2glytoucan = {}
    in_file = "downloads/glytoucan/current/export/matrixdb.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        matrixdb_id = "matrixdb:" + row[f_list.index("MatrixDBAccession")]
        matrixdbid2glytoucan[matrixdb_id] = glytoucan_ac


    log_file = path_obj["logs"] +  "%s_protein_matrixdb.log" % (species)
    FL = open(log_file, "w")
       

    newrow = ["uniprotkb_canonical_ac", "matrix_db_id", "matrix_db_label", "chebi_id", 
                "gag_interactor_biological_role", "protein_interactor_biological_role",
                "evidence", "saccharide","interaction_type", "source", "source_id", 
                "xref_key","xref_id", "src_xref_key","src_xref_id"] 
    print "\"%s\"" % ("\",\"".join(newrow))

    seen_row = {}
    data_frame = {}
    in_file = path_obj["downloads"] +  "/matrixdb/current/matrixdb_CORE.tab"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    for row in data_frame["data"]:
        ids_a = row[f_list.index("ID(s) interactor A")].replace("\"", "")
        ids_b = row[f_list.index("ID(s) interactor B")].replace("\"", "")
        alt_ids_a = row[f_list.index("Alt. ID(s) interactor A")].replace("\"", "")
        alt_ids_b = row[f_list.index("Alt. ID(s) interactor B")].replace("\"", "")
        pub_ids = row[f_list.index("Publication Identifier(s)")].replace("\"", "")
        roles_a = row[f_list.index("Biological role(s) interactor A")].replace("\"", "")
        roles_b = row[f_list.index("Biological role(s) interactor B")].replace("\"", "")
        interaction_type = "Glycosaminoglycan"

        ev_list = pub_ids.split("|") 
        if alt_ids_a.find("GAG") != -1 and ids_b.find("uniprotkb") != -1:
            matrixdb_id = alt_ids_a
            source_id = matrixdb_id
            source = "MatrixDB"
            ac = ids_b.split(":")[-1]
            chebi_id = ":".join(ids_a.split(":")[1:])
            gag_role = roles_a
            protein_role = roles_b
            if ac in ac2canon:
                canon = ac2canon[ac]
                for ev in ev_list:
                    saccharide = matrixdbid2glytoucan[matrixdb_id] if matrixdb_id in matrixdbid2glytoucan else ""
                    matrix_db_label = matrixdbid2label[matrixdb_id] if matrixdb_id in matrixdbid2label else ""
                    pmid = ev.split(":")[1]
                    matrixdb_id = matrixdb_id.replace("matrixdb:", "")
                    for pair in [["glycan_xref_matrixdb", matrixdb_id],["protein_xref_pubmed",pmid]]:

                        if pair[0] == "protein_xref_pubmed" and pair[1].isdigit() == False:
                            continue
                        newrow = [canon,matrixdb_id,matrix_db_label, chebi_id,gag_role,protein_role,ev,saccharide,interaction_type,source,source_id,pair[0], pair[1],"glycan_xref_matrixdb", matrixdb_id]
                        if saccharide == "" or saccharide in glycan_dict:
                            row_str = json.dumps(newrow)
                            if row_str not in seen_row:
                                seen_row[row_str] = True
                                print "\"%s\"" % ("\",\"".join(newrow))
                        else:
                            FL.write("%s,%s,glycan-not-in-masterlis\n" % (canon,saccharide))
        if alt_ids_b.find("GAG") != -1 and ids_a.find("uniprotkb") != -1:
            matrixdb_id = alt_ids_b
            source_id = matrixdb_id
            source = "PubMed"
            ac = ids_a.split(":")[-1]
            chebi_id = ":".join(ids_b.split(":")[1:])
            gag_role = roles_b
            protein_role = roles_a
            if ac in ac2canon:
                canon = ac2canon[ac]
                for ev in ev_list:
                    saccharide = matrixdbid2glytoucan[matrixdb_id] if matrixdb_id in matrixdbid2glytoucan else ""
                    matrix_db_label = matrixdbid2label[matrixdb_id] if matrixdb_id in matrixdbid2label else ""
                    pmid = ev.split(":")[1]
                    matrixdb_id = matrixdb_id.replace("matrixdb:", "")
                    for pair in [["glycan_xref_matrixdb", matrixdb_id],["protein_xref_pubmed",pmid]]:

                        if pair[0] == "protein_xref_pubmed" and pair[1].isdigit() == False:
                            continue
                        newrow = [canon,matrixdb_id,matrix_db_label,chebi_id,gag_role,protein_role,ev,saccharide,interaction_type,source,source_id,pair[0], pair[1], "glycan_xref_matrixdb", matrixdb_id]
                        if saccharide == "" or saccharide in glycan_dict:
                            row_str = json.dumps(newrow)
                            if row_str not in seen_row:
                                seen_row[row_str] = True
                                print "\"%s\"" % ("\",\"".join(newrow)) 
                        else:
                            FL.write("%s,%s,glycan-not-in-masterlis\n" %(canon,saccharide))              
    
    FL.close()

    return

def extract_glycogenes_glycoenzonto_ds(species):


    seen_row = {}
    sheet_obj = {}
    in_file = "downloads/glycogenes/current/%s_glycogenes_glycoenzonto.csv" % (species)
    libgly.load_sheet(sheet_obj, in_file, ",")
    f_list = sheet_obj["fields"]
    newrow = ["uniprotkb_canonical_ac", "gene_name", "gene_id", "recommended_name_full" ]
    print "\"%s\"" % ("\",\"".join(newrow))
    for row in sheet_obj["data"]:
        ac = row[f_list.index("UniProt")]
        gene_id = row[f_list.index("Entrez Gene ID")]
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        rec_name = canon2recname[canon] if canon in canon2recname else ""
        gene_name_list = canon2genenames[canon] if canon in canon2genenames else ""
        gene_id_list = canon2geneid[canon] if canon in canon2geneid else [""]
        for gene_name in gene_name_list:
            for gene_id in gene_id_list: 
                newrow = [canon, gene_name, gene_id, rec_name]
                row_s = json.dumps(newrow)
                if row_s not in seen_row:
                    print "\"%s\"" % ("\",\"".join(newrow))
                    seen_row[row_s] = True
 
    return 

def extract_genelocus_ds(species):
    
    data_grid = {"locusinfo":{}}
    sparqlutil.load_gene_locusinfo(data_grid, species)

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
            if gene_name in data_grid["locusinfo"]:
                ensg_id = data_grid["locusinfo"][gene_name]["ensgid"]
                chr_id = data_grid["locusinfo"][gene_name]["chrid"]
                start_pos = data_grid["locusinfo"][gene_name]["startpos"]
                end_pos = data_grid["locusinfo"][gene_name]["endpos"]
                strand = data_grid["locusinfo"][gene_name]["strand"]
                newrow = [canon, gene_name, ensg_id, chr_id, start_pos, end_pos, strand]
                print "\"%s\"" % ("\",\"".join(newrow))
        seen[canon] = True


    return


def extract_genelocus_ds_old(species):

    gtf_file = path_obj["downloads"] + "ucsc/gtf/%s.gtf" % (species)
    cmd = "readlink -f " + gtf_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")
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
   
    data_grid = {"ac2canon":{}, "interaction":{}}
    sparqlutil.load_interaction(data_grid, species)

    seen_row = {}
    newrow = ["uniprotkb_canonical_ac","intact_ac","participant_uniprotkb_ac","participant_intact_ac",
            "participant_uniprotkb_id","participant_gene_symbol", "participant_taxid","experiments"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for ac in data_grid["interaction"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            for o in data_grid["interaction"][ac]:
                newrow = [canon,o["intactac"], o["puniprotac"],o["pintactac"],o["puniprotid"], 
                        o["pgenename"], o["ptaxid"],o["experiments"]]
                row_str = json.dumps(newrow)
                if row_str not in seen_row:
                    seen_row[row_str] = True
                    print "\"%s\"" % ("\",\"".join(newrow))
    return



def extract_enzyme_annotation_uniprotkb_ds(species):

    data_grid = {"enzymeinfo":{}}
    sparqlutil.load_enzymelist_one(data_grid, species)
    sparqlutil.load_enzymelist_two(data_grid, species)

    seen_row = {}
    newrow = ["uniprotkb_canonical_ac", "enzyme_ec", "enzyme_activity"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for ac in data_grid["enzymeinfo"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            for o in data_grid["enzymeinfo"][ac]:
                newrow = [canon,o["ec"],o["activity"]]
                row_str = " ".join(newrow).lower()
                if row_str not in seen_row:
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen_row[row_str] = True


    return




def extract_glycoenzymes_ds(species, dataset):



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
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        if canon not in seen_canon:
            continue
        if len(row) != ncols:
            continue
        row[0] = canon
        print "\"%s\"" % ("\",\"".join(row))

    return





def extract_function_uniprotkb_ds(species):

    black_list = get_blacklisted_pmids(species)

    data_grid = {"function":{}}
    sparqlutil.load_function(data_grid, species)

    row = ["uniprotkb_canonical_ac","xref_key", "xref_id","annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    seen_row = {}
    for ac in data_grid["function"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            for ann_id in data_grid["function"][ac]:
                ann = data_grid["function"][ac][ann_id]["ann"]
                ann = ann.replace("\"", "`")
                xref_key = "protein_xref_uniprotkb_fun"
                row = [canon, xref_key,ac, ann]
                row_str = json.dumps(row)
                if row_str not in seen_row:
                    seen_row[row_str] = True
                    print "\"%s\""  % ("\",\"".join(row))
                for xref_id in data_grid["function"][ac][ann_id]["pmidlist"]:
                    xref_key = "protein_xref_pubmed"
                    if xref_id.isdigit() == False:
                        xref_key = "protein_xref_uniprotkb_fun"
                        if xref_id[0:3] == "MF_":
                            xref_key = "protein_xref_hamap"
                    if xref_id in black_list:
                        continue
                    row = [canon, xref_key, xref_id, ann]
                    row_str = json.dumps(row)
                    if row_str not in seen_row:
                        seen_row[row_str] = True
                        print "\"%s\""  % ("\",\"".join(row))

    
    return


def load_canon2geneid(species):

    in_file = path_obj["unreviewed"] + "%s_protein_xref_geneid.csv" % (species)
    if os.path.isfile(in_file) == False:
        return {}

    data_frame = {}
    libgly.load_sheet(data_frame,in_file, ",")
    f_list = data_frame["fields"]

    canon2geneid = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_id = row[f_list.index("xref_id")]
        if canon not in canon2geneid:
            canon2geneid[canon] = []
        if gene_id not in canon2geneid[canon] :
            canon2geneid[canon].append(gene_id)

    return canon2geneid


    
def extract_disease_glycosmos_ds(species):


    in_file = path_obj["downloads"] + "glycosmos/current/glycosmos_gdgdb.csv"
    data_frame = {}
    libgly.load_sheet(data_frame,in_file, ",")
    f_list = data_frame["fields"]
    geneid2ggdid = {}
    for row in data_frame["data"]:
        gene_id = row[f_list.index("gene_id")]
        ggd_id = row[f_list.index("dd_uri")].split("/")[-1].split("#")[-1]
        if gene_id not in geneid2ggdid:
            geneid2ggdid[gene_id] = []
        if ggd_id not in geneid2ggdid[gene_id]:
            geneid2ggdid[gene_id].append(ggd_id)


    in_file = path_obj["downloads"] + "glycosmos/current/graph-1/ggdonto.ttl"
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    ggdid2mimid = {}
    ggdid2name = {}
    with open(in_file, "r") as FR:
        for line in FR:
            if line[0:8] == "ggdonto:":
                ggd_id = line[8:].strip()
            if line.find("gdgsch:omimPhenoMIMnumber") != -1:
                mim_id = line.strip().split(" ")[1].replace("\"", "")
                if ggd_id not in ggdid2mimid:
                    ggdid2mimid[ggd_id] = []
                if mim_id not in ggdid2mimid[ggd_id]:
                    ggdid2mimid[ggd_id].append(mim_id)
            if line.find("ggdsch:diseasePreferredTerm") != -1:
                ggd_name =  " ".join(line.strip().split("@en")[0].split(" ")[1:]).replace("\"", "")
                ggdid2name[ggd_id] = ggd_name


    mondoid2doid = {}
    do_map = {}
    in_file = path_obj["unreviewed"] + "/protein_disease_idmap.csv"
    data_frame = {}
    libgly.load_sheet(data_frame,in_file, ",")
    f_list = data_frame["fields"]
    #['xref_id', 'xref_key', 'mondo_id', 'do_id', 'xref_map_src', 'do_map_src']
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        mondo_id = row[f_list.index("mondo_id")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key == "protein_xref_omim":
            mim_id = xref_id
            combo_id = "%s,%s" % (do_id, mondo_id)
            if mim_id not in do_map:
                do_map[mim_id] = []
            if combo_id not in do_map[mim_id]:
                do_map[mim_id].append(combo_id)
            

    row = ["uniprotkb_canonical_ac","xref_key", "xref_id","do_id", "mondo_id","mim_id",
            "ggd_id", "ggd_name"]
    print "\"%s\""  % ("\",\"".join(row))
    for canon in canon2geneid:
        for gene_id in canon2geneid[canon]:
            if gene_id not in geneid2ggdid:
                continue
            for ggd_id in geneid2ggdid[gene_id]:
                if ggd_id not in ggdid2mimid:
                    continue
                for mim_id in ggdid2mimid[ggd_id]:
                    ggd_name = ggdid2name[ggd_id] if ggd_id in ggdid2name else ""
                    if mim_id in do_map:
                        for combo_id in do_map[mim_id]:
                            do_id, mondo_id = combo_id.split(",")
                            row = [canon,"protein_xref_glycosmos_disease", gene_id,do_id,mondo_id,mim_id,ggd_id, ggd_name]
                            print "\"%s\""  % ("\",\"".join(row))
                    else:
                        do_id, mondo_id = "", ""
                        row = [canon,"protein_xref_glycosmos_disease", gene_id, do_id,mondo_id,
                                mim_id,ggd_id,ggd_name]
                        print "\"%s\""  % ("\",\"".join(row))


    return




def load_doid2name():

    tmp_dict = {}
    in_file = path_obj["unreviewed"] + "/protein_disease_names.csv"
    data_frame = {}
    libgly.load_sheet(data_frame,in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        name = row[f_list.index("name")]
        description = row[f_list.index("description")]
        if xref_key not in tmp_dict:
            tmp_dict[xref_key] = {}
        if xref_id not in tmp_dict[xref_key]:
            tmp_dict[xref_key][xref_id] = {"name":"", "description":"", "synonyms":[]}
        if row[f_list.index("name_type")] == "recommended_name":
            tmp_dict[xref_key][xref_id]["name"] = name
            tmp_dict[xref_key][xref_id]["description"] = description
        if row[f_list.index("name_type")] == "synonym" and name not in tmp_dict[xref_key][xref_id]["synonyms"]:
            tmp_dict[xref_key][xref_id]["synonyms"].append(name)
    

    return tmp_dict


def extract_disease_uniprotkb_ds(species):


    data_grid = {
        "ac2omim":{}, 
        "ac2mondo":{},
        "mimid2diseasename":{}
    }

    do_map, doid2mondoid, mondoid2doid, doid2mimid = {}, {}, {}, {}
    load_disease_id_map(do_map, doid2mondoid, mondoid2doid, doid2mimid)
           

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

    sparqlutil.load_mimid2disease_name(data_grid, species)
    sparqlutil.load_disease_mim(data_grid, species)


    seen = {}
    row = ["uniprotkb_canonical_ac","xref_key", "xref_id", "do_id", "mondo_id", "mim_id"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in ac2canon:
        canon = ac2canon[ac]
        if ac in data_grid["ac2mondo"]:
            for o in data_grid["ac2mondo"][ac]:
                combo_id = "%s|%s|%s" % ("", o["id"],"")
                if o["id"] in mondoid2doid:
                    combo_id = "%s|%s|%s" % (mondoid2doid[o["id"]], o["id"], "")
                row = [canon, "protein_xref_mondo", o["id"]] + combo_id.split("|")
                row_str = ",".join(row)
                if row_str not in seen:
                    print "\"%s\""  % ("\",\"".join(row))
                    seen[row_str] = True
        if ac in data_grid["ac2omim"]:
            for mim_id in data_grid["ac2omim"][ac]:
                d_name = ""
                if mim_id in data_grid["mimid2diseasename"]:
                    d_name = data_grid["mimid2diseasename"][mim_id]
                combo_id_list = ["|"]
                if mim_id in do_map["protein_xref_omim"]:
                    combo_id_list = do_map["protein_xref_omim"][mim_id]
                for combo_id in combo_id_list:
                    row = [canon, "protein_xref_omim", mim_id] + combo_id.split("|") + [mim_id]
                    row_str = ",".join(row)
                    if row_str not in seen:
                        print "\"%s\""  % ("\",\"".join(row))
                        seen[row_str] = True
    
    return


def load_disease_id_map(do_map, doid2mondoid, mondoid2doid, doid2mimid):

    in_file = path_obj["unreviewed"] + "/protein_disease_idmap.csv"
    data_frame = {}
    libgly.load_sheet(data_frame,in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        do_id = row[f_list.index("do_id")]
        mondo_id = row[f_list.index("mondo_id")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        if xref_key not in do_map:
            do_map[xref_key] = {}
        if xref_id not in do_map[xref_key]:
            do_map[xref_key][xref_id] = []
        combo_id = "%s|%s" % (do_id,mondo_id)
        if combo_id not in do_map[xref_key][xref_id]:
            do_map[xref_key][xref_id].append(combo_id)
        mondoid2doid[mondo_id] = do_id
        doid2mondoid[do_id] = mondo_id
        if do_id != "" and xref_key == "protein_xref_omim":
            mim_id = xref_id
            doid2mimid[do_id] = mim_id


    return


def extract_disease_genomics_england_ds(species):


    do_map, doid2mondoid, mondoid2doid, doid2mimid = {}, {}, {}, {}
    load_disease_id_map(do_map, doid2mondoid, mondoid2doid, doid2mimid)
           

    seen = {}
    row = ["uniprotkb_canonical_ac","xref_key", "xref_id", "do_id", "mondo_id", "mim_id"]
    print "\"%s\""  % ("\",\"".join(row))
               
    row_list = get_genomics_england_disease_rowlist(species)
    for row in row_list[1:]:
        do_id = row[-2]
        mondo_id = row[-1]
        mim_id, d_name = row[6].strip(), ""
        newrow = [row[0], "protein_xref_genomics_england", row[1],do_id,mondo_id,mim_id]
        newrow_str = ",".join(newrow)
        if newrow_str not in seen:
            print "\"%s\""  % ("\",\"".join(newrow))
            seen[newrow_str] = True
        continue

        #make row from mim_id
        combo_id_list = ["|"]
        if mim_id in do_map["protein_xref_omim"]:
            combo_id_list = do_map["protein_xref_omim"][mim_id]
        for combo_id in combo_id_list:
            newrow = [row[0], "protein_xref_omim", mim_id] + combo_id.split("|") + [mim_id]
            newrow_str = ",".join(newrow)
            if newrow_str not in seen:
                print "\"%s\""  % ("\",\"".join(newrow))
                seen[newrow_str] = True


def extract_disease_alliance_genome_ds(species):


    do_map, doid2mondoid, mondoid2doid, doid2mimid = {}, {}, {}, {}
    load_disease_id_map(do_map, doid2mondoid, mondoid2doid, doid2mimid)

    seen = {}
    row = ["uniprotkb_canonical_ac","xref_key", "xref_id", "do_id", "mondo_id", "mim_id"]
    print "\"%s\""  % ("\",\"".join(row))
    row_list = get_alliance_genome_disease_rowlist(species)
    for row in row_list:
        canon, do_id = row[0], row[1]
        mondo_id = doid2mondoid[do_id] if do_id in doid2mondoid else ""
        mim_id = doid2mimid[do_id] if do_id in doid2mimid else ""
        newrow = [canon, "protein_xref_genome_alliance", do_id,do_id,mondo_id, mim_id]
        newrow_str = ",".join(newrow)
        if newrow_str not in seen:
            print "\"%s\""  % ("\",\"".join(newrow))
            seen[newrow_str] = True

    return


def extract_disease_glyco_ds(species):

    return



def extract_disease_compiled_ds(dataset, species):


    data_frame = {}
    in_file = "compiled/%s_protein_%s.csv" % (species, dataset)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    print "\"%s\""  % ("\",\"".join(["uniprotkb_canonical_ac"] + f_list[1:]))
    log_file_one = path_obj["logs"] +   "%s_protein_%s.1.log" % (species, dataset)
    FL1 = open(log_file_one, "w")
    for row in data_frame["data"]:
        
        if dataset in ["disease_glygen"]:
            xref_key = "protein_xref_glygen_ds"
            ds = "%s_protein_%s" % (species, dataset)
            xref_id = ds2bco[ds] if ds in ds2bco else ""
            row[f_list.index("xref_key")] = xref_key
            row[f_list.index("xref_id")] = xref_id
        ac = row[0]
        if ac in ac2canon:
            canon = ac2canon[ac]
            newrow = [canon] + row[1:]
            print "\"%s\""  % ("\",\"".join(newrow))
        else:
            FL1.write("\"%s\"\n" % ("\",\"".join(row)))
    FL1.close()

    return


def extract_citations_generic_ds(species, dataset):


    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_%s.csv" % (species,dataset)
    in_file = in_file.replace("citations_", "")
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    
    biomarker_id_dict = {}
    seen_xref_id = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        src_xref_key, src_xref_id = "", ""
        if "src_xref_key" in f_list:
            src_xref_key = row[f_list.index("src_xref_key")]
            src_xref_id = row[f_list.index("src_xref_id")]
        if "dbsnp_id" in f_list:
            src_xref_key = "protein_xref_dbsnp"
            src_xref_id = row[f_list.index("dbsnp_id")]
        if "entrez_id" in f_list:
            src_xref_key = "protein_xref_entrez"
            src_xref_id = row[f_list.index("entrez_id")]
        if xref_key in ["protein_xref_pubmed", "protein_xref_doi"]:
            combo_id = "%s|%s|%s|%s|%s" % (canon, xref_key, xref_id, src_xref_key, src_xref_id)
            seen_xref_id[combo_id] = True
            if dataset == "citations_biomarkers":
                biomarker_id = row[f_list.index("biomarker_id")]
                if combo_id not in biomarker_id_dict:
                    biomarker_id_dict[combo_id] = {}
                biomarker_id_dict[combo_id][biomarker_id] = True


    newrow = ["uniprotkb_canonical_ac","title","journal_name",
            "publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    if dataset == "citations_biomarkers":
        newrow += ["biomarker_id"]
    print "\"%s\"" % ("\",\"".join(newrow))

    compiled_in_file = "compiled/doi_citations.csv"
    log_dict = {}
    log_file = path_obj["logs"] +  "%s_protein_%s.log" % (species,dataset)
    FL = open(log_file, "w")
    for combo_id in seen_xref_id:
        canon, xref_key, xref_id, src_xref_key, src_xref_id = combo_id.split("|")
        if src_xref_id == "":
            src_xref_key, src_xref_id = "protein_xref_pubmed", xref_id
        cite_info = {"row":[]}
        if xref_key == "protein_xref_pubmed":
            cite_info = libgly.get_citation(xref_id)
        elif xref_key == "protein_xref_doi":
            cite_info = libgly.get_doi_citation(xref_id, compiled_in_file)
        newrow = cite_info["row"]
        if newrow != []:
            if dataset == "citations_biomarkers":
                for biomarker_id in biomarker_id_dict[combo_id]:
                    print "\"%s\"" % ("\",\"".join([canon] + newrow + [xref_key, xref_id, src_xref_key, src_xref_id, biomarker_id]))
            else:
                print "\"%s\"" % ("\",\"".join([canon] + newrow + [xref_key, xref_id, src_xref_key, src_xref_id]))
        elif xref_id not in log_dict:
            FL.write("%s,%s\n" % (xref_id, ";".join(cite_info["flaglist"])))
            log_dict[xref_id] = True

    FL.close()



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
        pathway_id = row[f_list.index("xref_id")]
        if pathway_id not in pathway2canon:
            pathway2canon[pathway_id] = []
        if canon not in pathway2canon[pathway_id]:
            pathway2canon[pathway_id].append(canon)

    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_reactions_reactome.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    newrow = ["uniprotkb_canonical_ac","title","journal_name","publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    seen = {}
    log_dict = {}
    log_file = path_obj["logs"] +  "%s_protein_citations_reactome.log" % (species)
    FL = open(log_file, "w")
    

    for row in data_frame["data"]:
        reaction_id = row[f_list.index("reaction_id")].strip()
        pmid = row[f_list.index("pmid")].strip()
        pathway_id = row[f_list.index("pathway_id")].strip()
        if pmid in black_list or pmid == "":
            continue
        if pathway_id not in pathway2canon:
            if pmid not in log_dict:
                FL.write("%s,%s\n" % (pmid, "pathwayid='%s' not mapped to canonical" % (pathway_id)))
                log_dict[pmid] = True
            continue
        cite_info = libgly.get_citation(pmid)
        newrow = cite_info["row"]
        if newrow != []:
            for canon in pathway2canon[pathway_id]:
                combo_id = "%s %s" % (canon, pmid)
                if combo_id not in seen:
                    xref_key, xref_id = "protein_xref_pubmed", pmid
                    src_xref_key, src_xref_id = "protein_xref_reactome", pathway_id
                    out_row = [canon] + newrow + [xref_key, xref_id,src_xref_key, src_xref_id]
                    print "\"%s\"" % ("\",\"".join(out_row))
                    seen[combo_id] = True
        elif pmid not in log_dict:
            FL.write("%s,%s\n" % (pmid, ";".join(cite_info["flaglist"])))
            log_dict[pmid] = True

    FL.close()

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


    seen_combo = {}
    black_list = get_blacklisted_pmids(species)
    log_dict = {}
    log_file = path_obj["logs"] +  "%s_protein_citations_uniprotkb.log" % (species)
    FL = open(log_file, "w")
    
    row = ["uniprotkb_canonical_ac","title","journal_name", "publication_date", "authors"]
    row += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    print "\"%s\""  % ("\",\"".join(row))
   
    file_list = [path_obj["unreviewed"] +  "%s_protein_function_uniprotkb.csv" % (species)]
    file_list += [path_obj["unreviewed"] +  "%s_protein_ptm_annotation_uniprotkb.csv" % (species)] 
    file_list += [path_obj["unreviewed"] +  "%s_protein_site_annotation_uniprotkb.csv" % (species)]
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            xref_key,xref_id = row[f_list.index("xref_key")], row[f_list.index("xref_id")]
            if xref_key.find("xref_pubmed") != -1:
                pmid = xref_id
                if pmid.isdigit() == False:
                    continue
                if pmid in black_list:
                    continue
                combo = "%s|%s" % (canon, pmid)
                if combo in seen_combo:
                    continue
                seen_combo[combo] = True
                cite_info = libgly.get_citation(pmid)
                newrow = cite_info["row"]
                if newrow != []:
                    xref_key, xref_id = "protein_xref_pubmed", pmid
                    src_xref_key, src_xref_id = "protein_xref_uniprotkb", canon.split("-")[0]
                    print "\"%s\"" % ("\",\"".join([canon] + newrow + [xref_key, xref_id, src_xref_key, src_xref_id]))
                elif pmid not in log_dict:
                    FL.write("%s,%s\n" % (pmid, ";".join(cite_info["flaglist"])))
                    log_dict[pmid] = True

    data_grid = {"citelist":{}}
    sparqlutil.load_citelist(data_grid, species)
    for ac in data_grid["citelist"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            xref_key = "protein_xref_uniprotkb"
            xref_id = ac
            for obj in data_grid["citelist"][ac]:
                if obj["pmid"].isdigit() == False:
                    continue
                if obj["pmid"] in black_list:
                    continue
                combo = "%s|%s" % (canon, obj["pmid"])
                if combo in seen_combo:
                    continue
                seen_combo[combo] = True
                authors = ", ".join(obj["authorlist"])
                if authors == "":
                    continue
                xref_key, xref_id = "protein_xref_pubmed", obj["pmid"]
                src_xref_key, src_xref_id = "protein_xref_uniprotkb", ac
                row = [canon, obj["journaltitle"], obj["journalname"], 
                    obj["pubdate"], authors, xref_key, xref_id, src_xref_key, src_xref_id]
                print "\"%s\""  % ("\",\"".join(row))


    FL.close()
    return
 
 



def extract_canonicalsequences_ds(species):

    data_grid = {"ac2canon":{}, "genename":{}, "proteinnames":{},
            "proteinid":{}, "pe":{}, "isoformlist":{}, "isoforminfo":{}, 
            "isoformseq":{}, "isoformver":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)
    sparqlutil.load_isoformseq(data_grid, species)
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_protein_existence(data_grid, species)
    sparqlutil.load_protein_names(data_grid, species, "up:recommendedName")
    sparqlutil.load_protein_names(data_grid, species, "up:submittedName")
    sparqlutil.load_genename(data_grid, species)

    tax_name = species_obj[species]["long_name"]
    tax_id = species_obj[species]["tax_id"]

    seen_canon = {}
    for ac in data_grid["isoformlist"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            if canon in seen_canon:
                continue

            protein_id = data_grid["proteinid"][ac] if ac in data_grid["proteinid"] else ""
            pe = data_grid["pe"][ac] if ac in data_grid["pe"] else ""
            sv = data_grid["isoformver"][canon] if canon in data_grid["isoformver"] else ""
            gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
            full_name = ""
            if ac in data_grid["proteinnames"]:
                full_name = data_grid["proteinnames"][ac][0]["fullname"]

            seq = data_grid["isoformseq"][canon]
            status = data_grid["isoforminfo"][canon]["reviewed"]
            ac_lbl = "tr|%s|%s" % (canon,protein_id)
            isoform_id = canon.split("-")[-1]
            desc = "Isoform %s of %s OS=%s OX=%s GN=%s " %(isoform_id,full_name, tax_name, tax_id, gene_name)
            desc += "PE=%s " % (pe) if pe != "" else ""
            desc += "SV=%s " % (sv) if sv != "" else ""
            desc += "CANONICAL=%s" % (canon)

            if status == "1":
                ac_lbl = "sp|%s|%s" % (canon,protein_id)
            seq_obj = SeqRecord(Seq(seq,IUPAC.protein),id=ac_lbl, name=ac_lbl, description=desc)
            print "%s" % (seq_obj.format("fasta"))
            seen_canon[canon] = True




def extract_allsequences_ds(species):
    
    data_grid = {"ac2canon":{}, "genename":{}, "proteinnames":{}, 
            "proteinid":{}, "pe":{}, "isoformlist":{}, "isoforminfo":{}, 
            "isoformseq":{}, "isoformver":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)
    sparqlutil.load_isoformseq(data_grid, species)
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_protein_existence(data_grid, species)
    sparqlutil.load_protein_names(data_grid, species, "up:recommendedName")
    sparqlutil.load_protein_names(data_grid, species, "up:submittedName")
    sparqlutil.load_genename(data_grid, species)

    tax_name = species_obj[species]["long_name"]
    tax_id = species_obj[species]["tax_id"]
    seen_isoform = {}
    for ac in data_grid["isoformlist"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            protein_id = data_grid["proteinid"][ac] if ac in data_grid["proteinid"] else ""
            pe = data_grid["pe"][ac] if ac in data_grid["pe"] else ""
            gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
            full_name = ""
            if ac in data_grid["proteinnames"]:
                full_name = data_grid["proteinnames"][ac][0]["fullname"]
            for isoform in data_grid["isoformlist"][ac]:
                if isoform in seen_isoform:
                    continue

                seq = data_grid["isoformseq"][isoform]
                sv = data_grid["isoformver"][isoform] if isoform in data_grid["isoformver"] else ""
                status = data_grid["isoforminfo"][isoform]["reviewed"] 
                ac_lbl = "tr|%s|%s" % (isoform,protein_id)
                isoform_id = isoform.split("-")[-1]
                desc = "Isoform %s of %s OS=%s OX=%s GN=%s " %(isoform_id, full_name, tax_name, tax_id, gene_name)
                desc += "PE=%s " % (pe) if pe != "" else ""
                desc += "SV=%s " % (sv) if sv != "" else ""
                desc += "CANONICAL=%s" % (canon)                     

                if status == "1":
                    ac_lbl = "sp|%s|%s" % (isoform,protein_id)
                seq_obj = SeqRecord(Seq(seq,IUPAC.protein),id=ac_lbl, name=ac_lbl, description=desc)
                print "%s" % (seq_obj.format("fasta"))
                seen_isoform[isoform] = True                       



    return


def extract_sequenceinfo_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "genename":{}, "proteinnames":{},
            "proteinid":{}, "pe":{}, "isoformlist":{}, "isoforminfo":{},
            "isoformseq":{}, "isoformver":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)
    sparqlutil.load_isoformseq(data_grid, species)
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_protein_existence(data_grid, species)
    sparqlutil.load_protein_names(data_grid, species, "up:recommendedName")
    sparqlutil.load_protein_names(data_grid, species, "up:submittedName")
    sparqlutil.load_genename(data_grid, species)

    tax_name = species_obj[species]["long_name"]
    tax_id = species_obj[species]["tax_id"]
    
    newrow = ["uniprotkb_canonical_ac","uniprotkb_isoform_ac","protein_existence","sequence_version", "sequence_header"]
    print "\"%s\""  % ("\",\"".join(newrow))
    for ac in data_grid["isoformlist"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            protein_id = data_grid["proteinid"][ac] if ac in data_grid["proteinid"] else ""
            pe = data_grid["pe"][ac] if ac in data_grid["pe"] else ""
            gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
            full_name = ""
            if ac in data_grid["proteinnames"]:
                full_name = data_grid["proteinnames"][ac][0]["fullname"]
            for isoform in data_grid["isoformlist"][ac]:
                sv = data_grid["isoformver"][isoform] if isoform in data_grid["isoformver"] else ""
                status = data_grid["isoforminfo"][isoform]["reviewed"]
                desc = "tr|%s|%s " % (isoform,protein_id)
                if status == "1":
                    desc = "sp|%s|%s " % (isoform,protein_id)
                isoform_id = isoform.split("-")[-1]
                association_type = "Computationally mapped isoform of %s " % (full_name)
                if isoform.find(canon.split("-")[0]) != -1:
                    association_type = "Isoform %s of %s " % (isoform_id, full_name)
                desc += "%s OS=%s OX=%s GN=%s " % (association_type, tax_name, tax_id, gene_name)
                desc += "PE=%s " % (pe) if pe != "" else ""
                desc += "SV=%s " % (sv) if sv != "" else ""
                desc += "CANONICAL=%s" % (canon)

                newrow = [canon,isoform,str(pe),str(sv),desc]
                print "\"%s\""  % ("\",\"".join(newrow))            





def extract_info_uniprotkb_ds(species):

    data_grid = {"proteinid":{}, "isoformmass":{}, "isoformlen":{}}
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_isoformmass(data_grid, species)
    sparqlutil.load_isoformlen(data_grid, species)



    row = ["uniprotkb_canonical_ac","uniprotkb_id","uniprotkb_protein_mass","uniprotkb_protein_length"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in ac2canon:
        if ac == "":
            continue
        canon = ac2canon[ac]
        if canon != ac:
            continue
        ac = ac.split("-")[0]
        protein_id = data_grid["proteinid"][ac] if ac in data_grid["proteinid"] else ""
        canon_mass = data_grid["isoformmass"][canon] if canon in data_grid["isoformmass"] else -1
        canon_len = data_grid["isoformlen"][canon] if canon in data_grid["isoformlen"] else -1
        row = [canon, protein_id, canon_mass, canon_len]
        print "\"%s\""  % ("\",\"".join(row))



    return

def extract_ac2pdb_ds(species):

    ds_name = "ac2pdb"
    data_grid = {"ac2pdb":{}}
    sparqlutil.load_ac2pdb(data_grid, species)
    
    row = ["uniprotkb_canonical_ac","pdb_id"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in ac2canon:
        canon = ac2canon[ac]
        if ac in data_grid[ds_name]:
            for pdb_id in data_grid[ds_name][ac]:
                row = [canon, pdb_id]
                print "\"%s\""  % ("\",\"".join(row))


def extract_iedb_ds(species):


    fasta_file = "unreviewed/%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(fasta_file)

    f_map = {
        "iedb_id":"features.xrefs.id",
        "epitope_peptide":"features.epitopeSequence",
        "epitope_description":"features.description",
        "evidence_code":"features.evidences.code"
    }
    
    field_list_one = ["uniprotkb_canonical_ac"]
    field_list_two = ["iedb_id","epitope_peptide","epitope_description","evidence_code"]
    field_list_three = ["epitope_start_pos","epitope_end_pos","epitope_start_aa", "epitope_end_aa", "xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    print "\"%s\""  % ("\",\"".join(field_list_one + field_list_two + field_list_three)) 

    seen_row = {}
    tax_name = species_obj[species]["nt_file"].split("-proteome-")[1].replace(".nt", "")
    in_file = "downloads/ebi/current/IEDB-%s.tsv" % (tax_name)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("accession")]
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        if canon not in seq_hash:
            continue
        epitope_start_pos = int(row[f_list.index("features.begin")])
        epitope_end_pos = int(row[f_list.index("features.end")])
        epitope_peptide = row[f_list.index("features.epitopeSequence")]
        if epitope_start_pos < 0 or epitope_end_pos > len(seq_hash[canon]):
            continue
        peptide_in_canon = seq_hash[canon][epitope_start_pos-1:epitope_end_pos] 
        if epitope_peptide != peptide_in_canon:
            continue
 
        newrow = [canon]
        for f in field_list_two:
            val = row[f_list.index(f_map[f])]
            newrow.append(val)
        val_dict = {"xref_key":"protein_xref_pubmed", "src_xref_key":"protein_xref_iedb"}
        val_dict["epitope_start_pos"] = row[f_list.index("features.begin")]
        val_dict["epitope_end_pos"] = row[f_list.index("features.end")]
        val_dict["epitope_start_aa"] = peptide_in_canon[0]
        val_dict["epitope_end_aa"] = peptide_in_canon[-1]
        val_dict["xref_id"] = row[f_list.index("features.evidences.source.id")]
        val_dict["src_xref_id"] = row[f_list.index("features.xrefs.id")]
        for f in field_list_three:
            newrow.append(val_dict[f])
        row_str = json.dumps(newrow)
        if row_str in seen_row:
            continue
        seen_row[row_str] = True
        print "\"%s\""  % ("\",\"".join(newrow))
        
    return

def extract_xrefs_ds(species, ds_name):

    data_grid = {"reactome":{}}


    sheet_obj = {}
    if ds_name in ["xref_glycoprotdb"]:
        in_file = path_obj["downloads"] + "/glycoprotdb/protein_xref_glycoprotdb.csv"
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[f_list.index("uniprot_id")].split("-")[0]
            db_id = row[f_list.index("gpdb_id")]
            db_label = row[f_list.index("p_name")]
            o = {"id":db_id, "label":db_label}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            sheet_obj[ac].append(o)
    elif ds_name in ["xref_viralglycome"]:
        in_file = "compiled/%s_protein_xref_viralglycome.csv" % (species)
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[f_list.index("uniprot_canonical_ac")].split("-")[0]
            db_id = row[f_list.index("xref_id")]
            db_label = row[f_list.index("xref_label")]
            o = {"id":db_id, "label":db_label}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            sheet_obj[ac].append(o)
    elif ds_name in ["xref_iedb"]:
        tax_name = species_obj[species]["nt_file"].split("-proteome-")[1].replace(".nt", "")
        in_file = "downloads/ebi/current/IEDB-%s.tsv" % (tax_name)
        #print in_file
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, "\t")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac = row[f_list.index("accession")]
            db_id = row[f_list.index("features.xrefs.id")] 
            o = {"id":db_id, "label":"IEDB"}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            sheet_obj[ac].append(o)
    elif ds_name in ["xref_cfde_gdlpa"]:
        in_file = "compiled/%s_protein_xref_cfde_gdlpa.csv" % (species)
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            gene_name = row[f_list.index("gene_name")]
            canon = genename2canon[gene_name] if gene_name in genename2canon else ""
            ac = canon.split("-")[0]
            db_id = gene_name
            db_label = ""
            o = {"id":db_id, "label":db_label}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            sheet_obj[ac].append(o)
    elif ds_name in ["xref_pro_proteoform"]:
        sheet_obj = {}
        in_file = "downloads/pro/current/GlyGen_table_to_PRO.obo"
        cmd = "readlink -f " + in_file
        x = commands.getoutput(cmd)
        libgly.log_file_usage(x, "", "append")
        with open(in_file, "r") as FR:
            flag = False
            val_dict = {}
            for line in FR:
                line = line.strip()
                if line == "[Term]":
                    flag = True
                    if val_dict != {}:
                        ac = ""
                        if "is_a" in val_dict:
                            ac = val_dict["is_a"].strip().split(" ")[0].split(":")[1] 
                        db_id = val_dict["id"].strip()
                        if ac != "":
                            o = {"id":db_id, "label":"PRO"}
                            if ac not in sheet_obj:
                                sheet_obj[ac] = []
                            sheet_obj[ac].append(o)
                    val_dict = {}
                elif line == "":
                    flag = False
                if flag == True:
                    key = line.split(":")[0]
                    if key in ["id", "is_a"]:
                        val = ":".join(line.split(":")[1:])
                        val_dict[key] = val if key not in val_dict else val_dict[key] + "; " + val
            if val_dict != {}:
                ac = ""
                if "is_a" in val_dict:
                    ac = val_dict["is_a"].strip().split(" ")[0].split(":")[1]
                db_id = val_dict["id"].strip()
                if ac != "":
                    o = {"id":db_id, "label":"PRO"}
                    if ac not in sheet_obj:
                        sheet_obj[ac] = []
                    sheet_obj[ac].append(o)
    elif ds_name in ["xref_pharos"]:
        file_list = glob.glob(path_obj["downloads"] + "pharos/current/batch.*.json")
        for in_file in file_list:
            cmd = "readlink -f " + in_file
            x = commands.getoutput(cmd)
            libgly.log_file_usage(x, "", "append")
            doc = json.loads(open(in_file, "r").read())
            for doc in doc["content"]:
                db_id, db_label = doc["accession"], doc["name"]
                if db_id not in sheet_obj:
                    sheet_obj[db_id] = []
                o = {"id":db_id, "label":db_label}
                sheet_obj[db_id].append(o)
    elif ds_name in ["xref_uniprotkb"]:
        sparqlutil.log_nt_file_used(species)
        for ac in ac2canon:
            if ac.strip() == "":
                continue
            o = {"id":ac, "label":"UniProtKB"}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            sheet_obj[ac].append(o)
    elif ds_name in ["xref_unicarbkb"]:
        file_list = glob.glob(path_obj["unreviewed"] + "%s_proteoform_glycosylation_sites_*unicarbkb*.csv" % (species))
        seen_ref = {}
        for in_file in file_list:
            data_frame = {}
            libgly.load_sheet(data_frame, in_file, ",")
            f_list = data_frame["fields"]
            for row in data_frame["data"]:
                ac = row[f_list.index("uniprotkb_canonical_ac")].split("-")[0]
                src_xref_key = row[f_list.index("src_xref_key")]
                db_id = row[f_list.index("src_xref_id")]
                if db_id.strip() == "":
                    continue
                db_label = "UniCarbKB"
                o = {"ac":ac, "id":db_id, "label":db_label}
                if ac not in sheet_obj:
                    sheet_obj[ac] = []
                o_str = json.dumps(o)
                if o_str not in seen_ref:
                    sheet_obj[ac].append(o)
                    seen_ref[o_str] = True
    elif ds_name in ["xref_oglcnac_mcw", "xref_oglcnac_atlas", "xref_glyconnect"]:
        tmp_dict = {
            "xref_oglcnac_mcw":{"pref":"oglcnac_mcw", "lbl":"O-GlcNAc MCW"},
            "xref_oglcnac_atlas":{"pref":"oglcnac_atlas", "lbl":"O-GlcNAc Atlas"},
            "xref_glyconnect":{"pref":"glyconnect", "lbl":""}
        }
        pref = tmp_dict[ds_name]["pref"]
        db_label = tmp_dict[ds_name]["lbl"]
        in_file = path_obj["unreviewed"] + "%s_proteoform_glycosylation_sites_%s.csv" %(species,pref)
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        seen_ref = {}
        for row in data_frame["data"]:
            ac = row[f_list.index("uniprotkb_canonical_ac")].split("-")[0]
            db_id = row[f_list.index("src_xref_id")]
            if db_id.strip() == "":
                continue
            o = {"id":db_id, "label":db_label}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            o_str = json.dumps(o)
            if o_str not in seen_ref:
                sheet_obj[ac].append(o)
                seen_ref[o_str] = True
    elif ds_name == "xref_biomuta":
        in_file = path_obj["unreviewed"] + "%s_protein_mutation_cancer.csv" %(species)
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        seen_ref = {}
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            ac = canon.split("-")[0]
            o = {"id":ac, "label":""}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            o_str = json.dumps(o)
            if o_str not in seen_ref:
                sheet_obj[ac].append(o)
                seen_ref[o_str] = True
    elif ds_name == "xref_biomarkerkb":
        in_file = path_obj["unreviewed"] + "%s_protein_biomarkers.csv" %(species)
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        seen_ref = {}
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            ac = canon.split("-")[0]
            biomarker_id = row[f_list.index("biomarker_id")]
            o = {"id":biomarker_id, "label":""}
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            o_str = json.dumps(o)
            if o_str not in seen_ref:
                sheet_obj[ac].append(o)
                seen_ref[o_str] = True
    else:
        sparqlutil.load_ac2xref(sheet_obj, species, config_obj["xref"][ds_name])


    seen_row = {} 
    row = ["uniprotkb_canonical_ac","xref_key", "xref_id","xref_label"]
    print "\"%s\""  % ("\",\"".join(row))
    seen = {}
    xref_key = "protein_%s" % (ds_name)
    if ds_name == "xref_unicarbkb":
        xref_key = "protein_xref_unicarbkb_ds"
    for ac in ac2canon:
        canon = ac2canon[ac]
        if ac in sheet_obj:
            #consider one-to-one mapping only for xref_refseq
            if ds_name == "xref_refseq" and len(sheet_obj[ac]) > 1:
                selected_o = sheet_obj[ac][0]
                for o in sheet_obj[ac]:
                    if o["id"][0:3] == "NP_":
                        selected_o = o
                        break
                xref_id = selected_o["id"]
                xref_label = selected_o["label"].encode('ascii', 'ignore').decode('ascii')
                if ac in ac2canon_strict:
                    seen[canon] = True
                    row = [canon, xref_key, xref_id.split(".")[0], xref_label]
                    row_str = json.dumps(row)
                    if row_str not in seen_row:
                        print "\"%s\"" % ("\",\"".join(row))
                        seen_row[row_str] = True
            else:
                for o in sheet_obj[ac]:
                    xref_id = o["id"]
                    xref_label = o["label"].encode('ascii', 'ignore').decode('ascii')
                    if xref_key == "protein_xref_biomuta":
                        xref_id = ac
                    if ds_name == "xref_refseq":
                        if ac in ac2canon_strict:
                            row = [canon, xref_key, xref_id.split(".")[0], xref_label]
                            row_str = json.dumps(row)
                            if row_str not in seen_row:
                                print "\"%s\"" % ("\",\"".join(row))
                                seen_row[row_str] = True
                    elif ds_name == "xref_alphafolddb":
                        if ac in ac2canon_strict:
                            row = [canon, xref_key, xref_id, xref_label]
                            row_str = json.dumps(row)
                            if row_str not in seen_row:
                                print "\"%s\"" % ("\",\"".join(row))
                                seen_row[row_str] = True
                    else:
                        xref_id = ac if xref_key == "protein_xref_orthodb" else xref_id
                        row = [canon, xref_key, xref_id, xref_label]
                        row_str = json.dumps(row)
                        if row_str not in seen_row:
                            print "\"%s\"" % ("\",\"".join(row))
                            seen_row[row_str] = True
                    seen[canon] = True
        elif species.find("hcv") != -1 and ds_name in ["xref_hepatitisconline", "xref_viruspathogenresource"]:
            xref_id = "hcv"
            xref_label = ""
            row = [canon, xref_key, xref_id, xref_label]
            seen[canon] = True
            row_str = json.dumps(row)
            if row_str not in seen_row:
                print "\"%s\"" % ("\",\"".join(row))
                seen_row[row_str] = True

    #Add isoform level mappings for xref_refseq
    if ds_name == "xref_refseq":
        sheet_obj = {}
        sparqlutil.load_refseq2isoform(sheet_obj, species)
        for ac in ac2canon:
            canon = ac2canon[ac]
            if canon in seen:
                continue
            if canon in sheet_obj:
                selected = sheet_obj[canon][0]
                for refseq_isoform in sheet_obj[canon]:
                    if refseq_isoform[0:3] == "NP_":
                        selected = refseq_isoform
                        break
                xref_id = selected
                xref_label = ""
                row = [canon, xref_key, xref_id, xref_label]
                row_str = json.dumps(row)
                if row_str not in seen_row:
                    print "\"%s\"" % ("\",\"".join(row))
                    seen_row[row_str] = True
                if len(xref_id.split(".")) > 1:
                    row = [canon, xref_key, xref_id.split(".")[0], xref_label]
                    row_str = json.dumps(row)
                    if row_str not in seen_row:
                        print "\"%s\"" % ("\",\"".join(row))
                        seen_row[row_str] = True

    return



def extract_citations_refseq_ds(species):

    combo_id_list_one = get_refseq_pmidlist(species,"from_annotation")
    #combo_id_list_two = get_refseq_pmidlist(species,"from_references")
    combo_id_list_two = []
    combo_id_list = list(set(combo_id_list_one + combo_id_list_two))

    

    newrow = ["uniprotkb_canonical_ac","title","journal_name",
            "publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    print "\"%s\"" % ("\",\"".join(newrow))
    
    log_dict = {}
    log_file = path_obj["logs"] +  "%s_protein_citations_refseq.log" % (species)
    FL = open(log_file, "w")
    cite_info_dict = {}
    for combo_id in combo_id_list:
        canon, refseq_ac, pmid = combo_id.split(" ")
        xref_key = "protein_xref_refseq"
        xref_id = refseq_ac
        if pmid not in cite_info_dict:
            cite_info_dict[pmid] = libgly.get_citation(pmid)
        newrow = cite_info_dict[pmid]["row"]
        if newrow != []:
            xref_key, xref_id = "protein_xref_pubmed", pmid
            src_xref_key, src_xref_id = "protein_xref_refseq", refseq_ac
            print "\"%s\"" % ("\",\"".join([canon] + newrow + [xref_key, xref_id, src_xref_key, src_xref_id]))
        elif pmid not in log_dict:
            FL.write("%s,%s\n" % (pmid, ";".join(cite_info_dict[pmid]["flaglist"])))
            log_dict[pmid] = True

    FL.close()

    return




def get_refseq_pmidlist(species, target_src):
    
    black_list = get_blacklisted_pmids(species)

    seen = {}
    if target_src == "from_annotation":
        data_frame = {}
        in_file = path_obj["unreviewed"] +  "%s_protein_function_refseq.csv" % (species)
        if os.path.isfile(in_file) == False:
            return []
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            refseq_ac = row[f_list.index("refseq_ac")]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            if xref_key == "protein_xref_pubmed":
                pmid = xref_id
                #if pmid in black_list:
                #    continue
                combo_id = "%s %s %s" % (canon, refseq_ac, pmid)
                if combo_id not in seen:
                    seen[combo_id] = True
    elif target_src == "from_references":
        refseq2canon = {}
        canon2refseq = {}
        in_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
        if os.path.isfile(in_file) == False:
            return []
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            refseq = row[f_list.index("xref_id")]
            if canon not in canon2refseq:
                canon2refseq[canon] = refseq
                refseq2canon[refseq] = canon
        tax_id = species_obj[species]["tax_id"]
        in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_%s.gpff" % (tax_id)
        if species in ["sarscov1", "sarscov2", "hcv1", "hcv2"]:
            in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_viral.gpff"
        for record in SeqIO.parse(in_file, "genbank"):
            refseq_ac_one = record.id
            refseq_ac_two = record.id.split(".")[0]

            if refseq_ac_one not in refseq2canon and refseq_ac_two not in refseq2canon:
                continue
            refseq_ac = refseq_ac_one if refseq_ac_one in refseq2canon else refseq_ac_two
            canon = refseq2canon[refseq_ac]
            if "references" not in record.annotations:
                continue
            for ref in record.annotations["references"]:
                pmid = ref.pubmed_id
                if pmid != "":
                    #if pmid in black_list:
                    #    continue
                    combo_id = "%s %s %s" % (canon, refseq_ac, pmid)
                    if combo_id not in seen:
                        seen[combo_id] = True

    
    return seen.keys()

def load_glycan_masterdict():

    glycan_dict = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] +  "glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_ac")]
        glycan_dict[ac] = True

    return glycan_dict




def get_aa_dict():

    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_dict = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        one = row[f_list.index("one")]
        three = row[f_list.index("three")]
        aa_dict[three.upper()] = one.upper()
        aa_dict[one.upper()] = one.upper()

    return aa_dict



def get_glycosylation_sites(species):

    
    glycan_dict = load_glycan_masterdict()
    aa_dict = get_aa_dict()

    sites_dict = {}
    file_list = glob.glob("unreviewed/%s_proteoform_glycosylation_sites_*.csv" % (species))
    for in_file in file_list:
        if in_file.find("stat") != -1:
            continue
        source = in_file.split("_")[-1].split(".")[0]
        FR = open(in_file, "r") 
        idx = 0
        f_list = []
        for line in FR:
            idx += 1
            row = line.strip().split("\",\"")
            row[0] = row[0].replace("\"", "")
            row[-1] = row[-1].replace("\"", "")
            if idx == 1:
                f_list = row
                continue
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            gly_type = row[f_list.index("glycosylation_type")]
            continue_flag = False
            for f in ["glycosylation_site_uniprotkb", "amino_acid"]:
                if row[f_list.index(f)] == "":
                    continue_flag = True
            if continue_flag == True:
                continue
            pos = int(row[f_list.index("glycosylation_site_uniprotkb")].split("|")[0])
            amino_acid = row[f_list.index("amino_acid")].upper()
            ref_aa = aa_dict[amino_acid]
            xref_key = row[f_list.index("xref_key")]
            xref_id = row[f_list.index("xref_id")]
            glytoucan_ac = row[f_list.index("saccharide")].split("|")[0]
            score = 1
            if glytoucan_ac not in glycan_dict:
                glytoucan_ac = ""
                score = 0
            o = {"glytype":gly_type, "score":score, "xrefid":xref_id, "xrefkey":xref_key, 
                    "refaa":ref_aa, "glytoucanac":glytoucan_ac, "source":source}
            combo_id = "%s %s" % (canon, pos)
            if combo_id not in sites_dict:
                sites_dict[combo_id] = []
            sites_dict[combo_id].append(o)
        FR.close()

    return sites_dict


def get_mutation_sites(species):

    sites_dict = {}
    file_list = glob.glob("unreviewed/%s_protein_mutation_cancer.csv" % (species))
    for in_file in file_list:
        if in_file.find("stat") != -1:
            continue
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            #"uniprotkb_canonical_ac","aa_pos","ref_aa","alt_aa","chromosome_id","chr_pos","ref_nt","alt_nt","patients_positive","patients_tested","mut_freq","data_source","do_id","do_name","filter_flags"
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            pos = int(row[f_list.index("aa_pos")])
            ref_aa = row[f_list.index("ref_aa")]
            alt_aa = row[f_list.index("alt_aa")]
            do_id = row[f_list.index("do_id")]
            do_name = row[f_list.index("do_name")]
            flags = row[f_list.index("filter_flags")]
            score = flags.split(" ")[2].split("/")[0].replace("(", "")
            o = {"refaa":ref_aa, "altaa":alt_aa, "flags":flags, 
                    "score":score, "doid":do_id, "doname":do_name, "source":"biomuta"}
            combo_id = "%s %s" % (canon, pos)
            if combo_id not in sites_dict:
                sites_dict[combo_id] = []
            sites_dict[combo_id].append(o)


    return sites_dict


def get_mutation_effect_list(seq, aa_pos, ref_aa, alt_aa):

    e_list = []
    for offset in [0,1,2]:
        idx = aa_pos - offset
        if idx-1 < 0:
            continue
        ref_motif = seq[idx-1:idx+2]
        alt_motif = ""
        if offset == 0:
            alt_motif = alt_aa + ref_motif[offset+1:]
        else:
            alt_motif = ref_motif[:offset] + alt_aa + ref_motif[offset+1:]
        effect = "none"
        c_list = [is_n_glyco_motif(ref_motif)]
        c_list += [is_n_glyco_motif(alt_motif)]
        effect = "none"
        if c_list == [False,True]:
            e_list.append({"motifstart":idx,"refmotif":ref_motif, "altmotif":alt_motif,"effect":"n-glyco-sequon-gain"})
        elif c_list == [True,False]:
            e_list.append({"motifstart":idx,"refmotif":ref_motif, "altmotif":alt_motif,"effect":"n-glyco-sequon-loss"})
   
    if ref_aa in ["S","Y","T","K"] and alt_aa not in ["S","Y","T","K"]:
        e_list.append({"motifstart":aa_pos,"refmotif":ref_aa, "altmotif":alt_aa,"effect":"o-glyco-site-loss"})
    elif ref_aa in ["S","T"] and alt_aa not in ["S","T"]:
        e_list.append({"motifstart":aa_pos,"refmotif":ref_aa, "altmotif":alt_aa,"effect":"o-glyco-site-loss"})
    
    return e_list



def is_n_glyco_motif(motif):

    if len(motif) != 3:
        return False
    elif motif[0] == "N" and motif[1] != "P" and motif[2] in ["S", "T"]:
        return True
    else:
        return False



def check_motif_mutation_effect(ref_motif, alt_motif):

    #NXS/T(X!=P)
    if ref_motif[0] in ["N"]:
        if alt_motif[0] not in ["N"]:
            return "loss"
        elif alt_motif[1] in ["P"]:
            return "loss"
        elif alt_motif[2] not in ["S", "T"]:
            return "loss"
        else:
            return "none"
    else:
        return "N/A"

def extract_mutation_cancer_glycosylation_loss_ds(species):

    site_ann = {"gly":{}, "mut":{}}
    site_ann["gly"] = get_glycosylation_sites(species)
    site_ann["mut"] = get_mutation_sites(species)

    fasta_file = "unreviewed/%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(fasta_file)

    #set_one = set(site_ann["gly"].keys())
    #set_two = set(site_ann["mut"].keys())
    #site_list = list(set_one.intersection(set_two))
    newrow = [
        "uniprotkb_canonical_ac","mut_pos","mut_ref_aa","mut_alt_aa","ref_motif","alt_motif",
        "do_id","do_name","mut_flags_score","mut_flags","mut_source",
        "gly_pos","gly_ref_aa","saccharide","gly_xref_key","gly_xref_id","gly_flags_score","gly_source"
    ]
    print "\"%s\""  % ("\",\"".join(newrow))

    seen_row = {}
    site_list = site_ann["gly"].keys()
    for combo_id in site_list:
        canon = combo_id.split(" ")[0]
        gly_pos = int(combo_id.split(" ")[1])
        g_obj_list = []
        m_obj_list = []
        for g_obj in site_ann["gly"][combo_id]:
            #if g_obj["glytoucanac"] == "":
            #    continue
            if g_obj["refaa"] not in ["N"]:
                continue
            gly_ref_aa = g_obj["refaa"]
            for offset in [0,1,2]:
                mut_pos = gly_pos + offset
                idx = gly_pos - 1
                ref_motif = seq_hash[canon][idx:idx+3]
                offset_combo_id = "%s %s" % (canon, mut_pos)
                if offset_combo_id not in site_ann["mut"]:
                    continue
                for m_obj in site_ann["mut"][offset_combo_id]:
                    mut_ref_aa, mut_alt_aa = m_obj["refaa"],m_obj["altaa"]
                    alt_motif = ""
                    if offset == 0:
                        alt_motif = m_obj["altaa"] + ref_motif[offset+1:]
                    else:
                        alt_motif = ref_motif[:offset] + m_obj["altaa"] + ref_motif[offset+1:]
                    effect = check_motif_mutation_effect(ref_motif, alt_motif)
                    if effect in ["loss"]:
                        newrow = [
                            canon,
                            str(mut_pos),mut_ref_aa,mut_alt_aa,ref_motif,alt_motif,
                            m_obj["doid"],m_obj["doname"],str(m_obj["score"]),
                            m_obj["flags"], m_obj["source"],
                            str(gly_pos),gly_ref_aa,g_obj["glytoucanac"],
                            g_obj["xrefkey"],g_obj["xrefid"],str(g_obj["score"]),
                            g_obj["source"]
                        ]
                        row_str = json.dumps(newrow)
                        if row_str not in seen_row:
                            seen_row[row_str] = True
                            print "\"%s\""  % ("\",\"".join(newrow))
                       

    return



def extract_integrated_site_annotation_ds(species):
   
    data_grid = {"genename":{},
            "proteinnames":{},"shortname":{}, "mutagenann":{}}

    sparqlutil.load_genename(data_grid, species)
    sparqlutil.load_mutagen_annotation(data_grid, species)
    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_canonicalsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())
    site_ann = {"gly":{}, "pho":{}, "mutagen":{}, "dbsnp":{}, "ptm":{}}
    for ac in data_grid["mutagenann"]:
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        for o in data_grid["mutagenann"][ac]:
            start_pos = o["startpos"]
            end_pos = o["endpos"]
            if start_pos != end_pos:
                continue
            ecoid = o["ecoid"]
            pmid = o["pmid"]
            if pmid == "":
                continue
            for pos in xrange(start_pos, end_pos+1):
                combo_id = "%s %s" % (canon, pos)
                if combo_id not in site_ann["mutagen"]:
                    site_ann["mutagen"][combo_id] = []
                site_ann["mutagen"][combo_id].append(o)

            
    site2glytoucan = {}
    site2glytype = {}
    file_list = glob.glob(path_obj["unreviewed"] + "%s_proteoform_glycosylation_sites_*.csv" % (species))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            gly_type = row[f_list.index("glycosylation_type")]
            pos = int(row[f_list.index("glycosylation_site_uniprotkb")])
            pmid = row[f_list.index("evidence")]
            if pmid == "":
                continue
            glytoucan_ac = row[f_list.index("saccharide")]
            o = {"glytype":gly_type, "pmid":pmid, "glytoucanac":glytoucan_ac}
            combo_id = "%s %s" % (canon, pos)
            if glytoucan_ac.strip() != "":
                site2glytoucan[combo_id] = glytoucan_ac
            if combo_id not in site2glytype:
                site2glytype[combo_id] = []
            if gly_type.lower() not in site2glytype[combo_id]:
                site2glytype[combo_id].append(gly_type.lower())
            if combo_id not in site_ann["gly"]:
                site_ann["gly"][combo_id] = []
            site_ann["gly"][combo_id].append(o)

    file_list = glob.glob(path_obj["unreviewed"] + "%s_proteoform_phosphorylation_sites_*.csv" % (species))
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            amino_acid = row[f_list.index("amino_acid")]
            residue = row[f_list.index("phosphorylated_residue")]
            pos = int(row[f_list.index("phosphorylation_site_uniprotkb")])
            pmid = row[f_list.index("evidence")]
            ecoid = row[f_list.index("eco_id")]
            if pmid == "":
                continue
            o = {"aminoacid":amino_acid, "residue":residue,  
                        "pmid":pmid, "eco_id":ecoid}
            combo_id = "%s %s" % (canon, pos)
            if combo_id not in site_ann["pho"]:
                site_ann["pho"][combo_id] = []
            site_ann["pho"][combo_id].append(o)

    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_three_list = []
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        aa_three_list.append(three)

    kw_list = ["glycosylated", "glycosylation", "phosphorylated","phosphorylation"]
    in_file = path_obj["unreviewed"]+"%s_protein_ptm_annotation_uniprotkb.csv"%(species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        pmid = row[f_list.index("pmid")]
        ecoid = row[f_list.index("eco_id")]
        ann = row[f_list.index("ptm_annotation")]

        pos_dict = {}
        word_list = ann.lower().split(" ")
        for j in xrange(0, len(word_list)):
            if word_list[j] in kw_list:
                if j + 2 < len(word_list):
                    if word_list[j+2].find("-") == -1:
                        continue
                    aa = word_list[j+2].split("-")[0]
                    aa = aa[0].upper() + aa[1:]
                    if aa in aa_three_list:
                        for badchar in [".", ",", ";",":"]:
                            word_list[j+2] = word_list[j+2].replace(badchar, "")
                        pos_str_list = word_list[j+2].split("-")[1].split("/")
                        for pos_str in pos_str_list:
                            if pos_str.isdigit() == False:
                                continue
                            pos = int(pos_str)
                            kw = word_list[j]
                            if kw not in pos_dict:
                                pos_dict[kw] = []
                            pos_aa_combo = "%s %s" % (pos, aa)
                            if pos_aa_combo not in pos_dict[kw]:
                                pos_dict[kw].append(pos_aa_combo)
        for kw in pos_dict:
            for pos_aa_combo in pos_dict[kw]:
                pos, aa = pos_aa_combo.split(" ")[0], pos_aa_combo.split(" ")[1] 
                combo_id = "%s %s" % (canon, pos)
                kw_new = kw
                kw_new = "Phosphorylation" if "phosphorylat" in kw else kw_new
                kw_new = "Glycosylation" if "glycosylat" in kw else kw_new
                o = {"kw":kw, "ecoid":ecoid, "aminoacid":aa, "comment":ann}
                if combo_id not in site_ann["ptm"]:
                    site_ann["ptm"][combo_id] = []
                site_ann["ptm"][combo_id].append(o)


    tax_name = species_obj[species]["long_name"].lower().replace(" ", "-")
    in_file = path_obj["downloads"] + "/ebi/current/dbSNP-%s.tsv" % (tax_name)
    #in_file = "tmp/dbSNP-homo-sapiens.tsv"

    f_list = []
    k_list = ["uniprotkb_accession", "data_source", "dbsnp_id", "cosmic_id",
        "ref_aa", "alt_aa", "mutation_type", "polyphen_prediction","somatic_status",
        "sift_prediction", "disease", "begin_aa_pos", "end_aa_pos"]
    
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    import sys
    reload(sys)
    sys.setdefaultencoding('utf-8')
    import io
    with io.open(in_file, "r", encoding="utf-8",errors="ignore") as FR:
        lcount = 0
        for line in FR:
            lcount += 1
            #if lcount%1000000 == 0:
            #    print lcount
            row = line.split("\t")
            if lcount == 1:
                f_list = row
            else:
                o = {}
                for k in k_list:
                    o[k] = row[f_list.index(k)]
                ac = row[0].split("-")[0]
                if ac not in ac2canon:
                    continue
                canon = ac2canon[ac]
                start_pos = int(o["begin_aa_pos"])
                for offset in xrange(0, 3):
                    combo_id = "%s %s" % (canon, start_pos - offset)
                    for k_one in ["gly", "mutagen", "pho"]:
                        if combo_id in site_ann[k_one]:
                            if combo_id not in site_ann["dbsnp"]:
                                site_ann["dbsnp"][combo_id] = []
                            o["offset"] = offset
                            site_ann["dbsnp"][combo_id].append(o)
                            break


    id_list = []
    for k in site_ann:
        id_list += site_ann[k].keys()


    newrow = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
                "combination_flag", "gly_type", "has_glytoucan",
            "annotation_type","annotation"
    ]
    print "\"%s\""  % ("\",\"".join(newrow))

    seen_row = {}
    for combo_id in list(set(id_list)):
        canon, pos = combo_id.split(" ")
        flag_list = []
        row_list = []
        for k_one in site_ann.keys():
            if combo_id in site_ann[k_one]:
                flag_list.append(k_one)
                for o in site_ann[k_one][combo_id]:
                    ann_list = []
                    for k_two in o:
                        ann_list.append("%s=%s" % (k_two, o[k_two]))
                    ann = ";".join(ann_list).replace("\"", "")
                    row_list.append([canon, str(pos), k_one, ann])

        if len(flag_list) > 1:
            for row in row_list:
                combo_flag = ";".join(sorted(flag_list))
                has_glytoucan = "yes" if combo_id in site2glytoucan else "no"
                gly_type = "none"
                if combo_id in site2glytype:
                    gly_type = ";".join(sorted(site2glytype[combo_id]))
                newrow = row[0:2] + [combo_flag, gly_type, has_glytoucan] + row[2:]
                xref_key, xref_id = newrow[-2], newrow[-1]
                if xref_key == "protein_xref_pubmed" and xref_id.isdigit() == False:
                    continue 
                row_str = " ".join(newrow).lower()
                if row_str not in seen_row:
                    print "\"%s\""  % ("\",\"".join(newrow))
                    seen_row[row_str] = True


    return








def extract_signalp_annotation_ds(species):
   
    data_grid = {"genename":{},
                "proteinnames":{},"shortname":{}, "signalpann":{}}
    sparqlutil.load_protein_names(data_grid, species, "up:recommendedName")
    sparqlutil.load_genename(data_grid, species)
    sparqlutil.load_signalp_annotation(data_grid, species)

    seq_hash = {}
    fasta_file = path_obj["unreviewed"] +  "%s_protein_canonicalsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())

    row = ["uniprotkb_canonical_ac","gene_symbol","uniprotkb_protein_recommended_name",
        "uniprotkb_protein_length","uniprotkb_protein_signalp_length",
        "uniprotkb_protein_signalp_cleaved_length","signalp_start_pos",
        "signalp_end_pos","uniprotkb_protein_sequence","uniprotkb_protein_signalp_sequence",
        "uniprotkb_protein_signalp_cleaved_sequence","eco_id","evidence"]
    print "\"%s\""  % ("\",\"".join(row))

    seen_row = {}    
    
    for ac in data_grid["signalpann"]:
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        full_seq = seq_hash[canon]
        protein_len = str(len(full_seq))
        full_name = ""
        if ac in data_grid["proteinnames"]:
            obj = data_grid["proteinnames"][ac][0]
            full_name = obj["fullname"] if "fullname" in obj else ""


        for o in data_grid["signalpann"][ac]:
            start_pos = o["startpos"]
            end_pos = o["endpos"]
            ecoid = o["ecoid"]
            pmid = o["pmid"]
            peptide_seq = full_seq[start_pos-1:end_pos]
            cleaved_seq = full_seq[end_pos:]
            row = [canon, gene_name, full_name,str(len(full_seq)),
                    str(len(peptide_seq)), str(len(cleaved_seq)),
                    str(start_pos),str(end_pos),full_seq,peptide_seq,cleaved_seq,ecoid,pmid]
            row_str = json.dumps(row)
            if row_str not in seen_row:
                seen_row[row_str] = True
                print "\"%s\""  % ("\",\"".join(row))


    return

def extract_signalp_sequences_ds(species, seq_type):

    data_grid = {"genename":{},
                "proteinnames":{},"shortname":{}, "signalpann":{}}
    sparqlutil.load_protein_names(data_grid, species, "up:recommendedName")
    sparqlutil.load_genename(data_grid, species)
    sparqlutil.load_signalp_annotation(data_grid, species)

    seq_hash = {}
    fasta_file = path_obj["unreviewed"] +  "%s_protein_canonicalsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())

    tax_name = species_obj[species]["long_name"]
    tax_id = species_obj[species]["tax_id"]

    for ac in data_grid["signalpann"]:
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        full_seq = seq_hash[canon]
        full_name = ""
        if ac in data_grid["proteinnames"]:
            obj = data_grid["proteinnames"][ac][0]
            full_name = obj["fullname"] if "fullname" in obj else ""
        for o in data_grid["signalpann"][ac]:
            start_pos = o["startpos"]
            end_pos = o["endpos"]
            peptide_seq = full_seq[start_pos-1:end_pos]
            cleaved_seq = full_seq[end_pos:]
            ac_lbl = "%s" % (canon)
            if seq_type == "peptide_seq":
                desc = "%s OS=%s OX=%s GN=%s (signal peptide %s-%s)" %(full_name, tax_name, tax_id, gene_name, start_pos, end_pos)
                seq = peptide_seq
            elif seq_type == "cleaved_seq":
                desc = "%s OS=%s OX=%s GN=%s (cleaved signalp sequence %s-%s)" %(full_name, tax_name, tax_id, gene_name, end_pos + 1, len(full_seq))
                seq = cleaved_seq
            elif seq_type == "full_seq":
                desc = "%s OS=%s OX=%s GN=%s (full signalp sequence)" %(full_name, tax_name, tax_id, gene_name)
                seq = cleaved_seq
            seq_obj = SeqRecord(Seq(seq,IUPAC.protein),id=ac_lbl,name=ac_lbl,description=desc)
            print "%s" % (seq_obj.format("fasta"))

    return


def extract_proteinnames_ds(species, ds_name):


    data_grid = {"proteinnames":{}, "shortname":{}}

    
    predicate = "up:recommendedName"
    row = ["uniprotkb_canonical_ac","recommended_name_full","recommended_name_short", "ec_name"]
    if ds_name == "altnames":
        row = ["uniprotkb_canonical_ac","alternative_name_full","alternative_name_short", 
                "ec_name"]
        predicate = "up:alternativeName"
    elif ds_name == "submittednames":
        row = ["uniprotkb_canonical_ac","submitted_name_full","submitted_name_short", 
                "ec_name"]
        predicate = "up:submittedName"

    sparqlutil.load_protein_names(data_grid, species, predicate)

    print "\"%s\""  % ("\",\"".join(row))
    for ac in ac2canon:
        canon = ac2canon[ac]
        if ac != canon.split("-")[0]:
            continue
        full_name, short_name = "", ""
        if ac in data_grid["proteinnames"]:
            for o in data_grid["proteinnames"][ac]:
                row = [canon]
                for k in ["fullname", "shortname", "ecname"]:
                    row.append(o[k])
                print "\"%s\""  % ("\",\"".join(row))

    return

def extract_genenames_uniprotkb_ds(species):


    seen_row = {}

    data_grid = {"genenames":{}}
    sparqlutil.load_gene_names(data_grid, species)
    row = ["uniprotkb_canonical_ac", "gene_symbol_recommended", "gene_symbol_alternative", "orf_name"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["genenames"]:
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        for o in data_grid["genenames"][ac]:
            if list(set(o.values())) == [""]:
                continue
            row = [canon]
            for k in ["preflabel", "altlabel", "orfname"]:
                row.append(o[k])
            row_str = json.dumps(row)
            if row_str not in seen_row:
                seen_row[row_str] = True
                print "\"%s\""  % ("\",\"".join(row))


def extract_pdb_shortlist_ds(species):

    data_grid = {"pdbinfo":{}}
    sparqlutil.load_pdbinfo(data_grid, species)
    row = ["uniprotkb_canonical_ac","pdb_id", "method", "resolution"]
    print "\"%s\""  % ("\",\"".join(row))

    seen = {}
    selected_dict = {"xray":{}, "nmr":{}}
    for ac in ac2canon:
        canon = ac2canon[ac]
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

    data_grid = {"ptmann":{}}
    sparqlutil.load_ptm_annotation(data_grid, species)
    row = ["uniprotkb_canonical_ac","xref_key","xref_id", "eco_id", "ptm_annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    seen = {}
    for ac in ac2canon:
        canon = ac2canon[ac]
        if ac in data_grid["ptmann"]:
            for obj in data_grid["ptmann"][ac]:
                if obj["pmid"].strip() == "" or obj["pmid"].isdigit() == False:
                    continue
                for row in [
                    [canon,"protein_xref_uniprotkb", ac, obj["ecoid"], obj["comment"]]
                    ,[canon,"protein_xref_pubmed", obj["pmid"], obj["ecoid"], obj["comment"]]
                ]:
                    row_str = "\",\"".join(row)
                    if row_str not in seen:
                        print "\"%s\""  % (row_str)
                        seen[row_str] = True

    return


def extract_ntdata_ds(species):
    nt_file = species_obj[species]["nt_file"]
    ebi_dir = "downloads/ebi/current/"
    #ebi_dir = "downloads/ebi/2020_01/"

    in_file = "%s%s" % (ebi_dir,nt_file)
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")


    cmd = "cp %s unreviewed/%s_protein_ntdata.nt" % (in_file, species)
    x = commands.getoutput(cmd)
    return



def extract_site_annotation_uniprotkb_ds(species):


    black_list = get_blacklisted_pmids(species)

    canon2genename = load_canon2genename(species)

    data_grid = {"genename":{}, "siteann":{}}
    sparqlutil.load_genename(data_grid, species)
    
    seq_hash = {}
    fasta_file = path_obj["unreviewed"] + "%s_protein_canonicalsequences.fasta" % (species)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())


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
        "Signal_Peptide_Annotation",
        "Mutagenesis_Annotation", 
        "Natural_Variant_Annotation"
    ]
    for ann_type in ann_type_list:
        sparqlutil.load_site_annotation(data_grid, species, ann_type)
   
    row = ["uniprotkb_canonical_ac","gene_symbol", "ann_type", "start_pos","end_pos", "annotation", "ecoid","ref_aa","alt_aa","xref_key","xref_id"]
    print "\"%s\""  % ("\",\"".join(row))


    seen_row = {} 
    for ac in ac2canon:
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        gene_name = canon2genename[canon] if canon in canon2genename else ""
        if ac in data_grid["siteann"]:
            for obj in data_grid["siteann"][ac]:
                pmid = obj["pmid"] if "pmid" in obj else ""
                ecoid = obj["ecoid"] if "ecoid" in obj else ""
                ref_aa, alt_aa = "", ""
                if len(seq_hash[canon]) >= int(obj["end"]):
                    ref_aa = seq_hash[canon][int(obj["start"])-1:int(obj["end"])]
                    alt_aa = obj["subs"] if "subs" in obj else ""
                if ref_aa == "":
                    continue

                row = [canon, gene_name,obj["anntype"], str(obj["start"]), str(obj["end"]),
                        obj["ann"], ecoid,ref_aa,alt_aa,"protein_xref_uniprotkb",ac]
                row_str = json.dumps(row)
                if row_str not in seen_row:
                    print "\"%s\""  % ("\",\"".join(row))
                    seen_row[row_str] = True

                if pmid != "" and pmid.isdigit():
                    if pmid in black_list:
                        continue
                    row = [canon, gene_name,obj["anntype"], str(obj["start"]), str(obj["end"]),
                        obj["ann"], ecoid,ref_aa,alt_aa,"protein_xref_pubmed",pmid]
                    row_str = json.dumps(row)

                    if row_str not in seen_row:
                        print "\"%s\""  % ("\",\"".join(row))
                        seen_row[row_str] = True

    return




def extract_pro_annotation_ds(species):


    in_file = "downloads/pro/current/GlyGen_table_to_PRO.obo"
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    obj_list = []
    with open(in_file, "r") as FR:
        flag = False
        val_dict = {}
        for line in FR:
            line = line.strip()
            if line == "[Term]":
                flag = True
                if val_dict != {}:
                    obj_list.append(val_dict)
                val_dict = {}
            elif line == "":
                flag = False
            if flag == True:
                key = line.split(":")[0]
                val = ":".join(line.split(":")[1:])
                val_dict[key] = val if key not in val_dict else val_dict[key] + "; " + val
        if val_dict != {}:
            obj_list.append(val_dict)
    
    field_dict = {
        "id":"pro_id",
        "name":"pro_protein_name",
        "def":"pro_protein_definition",
        "comment":"pro_comment",
        "synonym":"pro_synonym",
        "relationship":"pro_relationship"
    }
    f_list = ["id", "name", "def", "comment", "synonym", "relationship"]
    newrow = ["uniprotkb_canonical_ac"]
    for f in f_list:
        newrow.append(field_dict[f])
    newrow += ["xref_key", "xref_id"]
    print "\"%s\""  % ("\",\"".join(newrow))
    for obj in obj_list:
        newrow = []
        for f in f_list:
            val = obj[f].strip().replace("\"", "") if f in obj else ""
            newrow.append(val)
        if newrow[0] == "" or newrow[2] == "":
            continue
        ac = obj["is_a"].strip().split(" ")[0].split(":")[1] if "is_a" in obj else ""
        if ac not in ac2canon:
            continue
        newrow = [ac2canon[ac]] + newrow + ["protein_xref_pro_proteoform", obj["id"].strip()]
        print "\"%s\""  % ("\",\"".join(newrow))

    return



def extract_go_annotation_ds(species):

    data_grid = {"genename":{}, "goann":{}}
    sparqlutil.load_go_annotation(data_grid, species)
    sparqlutil.load_genename(data_grid, species)

    seen_row = {}
    row = ["uniprotkb_canonical_ac","gene_symbol", "go_term_id","go_term_label", "go_term_category", "eco_id", "pmid"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in ac2canon:
        canon = ac2canon[ac]
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        if ac in data_grid["goann"]:
            for obj in data_grid["goann"][ac]:
                row = [canon, gene_name, obj["goid"], obj["goterm"],obj["gocat"],
                        obj["ecoid"],obj["pmid"]]
                row_str = json.dumps(row)
                if row_str not in seen_row:
                    seen_row[row_str] = True
                    print "\"%s\""  % ("\",\"".join(row))

    return





def extract_transcriptlocus_ds(species):


    data_grid = {"ac2canon":{}, "isoformlist":{}, "isoforminfo":{}, "locusinfo":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_transcript_locusinfo(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)


    seen = {}
    row = ["uniprotkb_canonical_ac","uniprotkb_isoform_ac","transcript_id","peptide_id","chromosome_id",
            "start_pos","end_pos","strand"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["isoformlist"]:
        if ac in ac2canon:
            canon = ac2canon[ac]
            for isoform in data_grid["isoformlist"][ac]:
                #print "Robel", canon, isoform, isoform in data_grid["locusinfo"]
                if isoform in data_grid["locusinfo"]:
                    o = data_grid["locusinfo"][isoform]
                    row = [canon,isoform,o["trsid"], o["pepid"],
                            str(o["chrid"]),str(o["startpos"]),str(o["endpos"]),str(o["strand"])]
                    row_str = json.dumps(row)
                    if row_str not in seen:
                        print "\"%s\"" % ("\",\"".join(row))
                        seen[row_str] = True



def extract_masterlist_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "genename":{}, "isoforminfo":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_genename(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)

    row = ["uniprotkb_canonical_ac","status","gene_name","reviewed_isoforms","unreviewed_isoforms"]
    print "\"%s\""  % ("\",\"".join(row))

    for ac in data_grid["isoformlist"]:
        #print "Robel-1", ac, ac in data_grid["ac2canon"]
        gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
        list_one, list_two = [], []
        for isoform in data_grid["isoformlist"][ac]:
            #print "Robel-2", ac, isoform
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
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        refseq = row[f_list.index("xref_id")]
        canon2refseq[canon] = refseq
        refseq2canon[refseq] = canon

    data_frame = {}
    field = "xxx"
    in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_%s.gpff" % (tax_id)
    if species in ["sarscov1", "sarscov2", "hcv1", "hcv2"]:
        in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_viral.gpff"

    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    with open(in_file, "r") as FR:
        for line in FR:
            newfield = line[0:12].strip() 
            field = newfield if len(newfield) > 0 else field
            value = line[12:].strip()
            if field == "ACCESSION":
            #if field == "VERSION":
                tmp_ac_list = line.replace("ACCESSION", "").strip().split(" ")
                for tmp_ac in tmp_ac_list:
                    ac = tmp_ac
                    if ac.find(".") != -1:
                        ac = ".".join(tmp_ac.strip().split(".")[:-1])
                    if ac in refseq2canon:
                        data_frame[ac] = {"ac":ac, "summary":""}
                        flag = False
                        break
            elif field == "COMMENT":
                if ac in refseq2canon:
                    if line.strip()[0:8] == "Summary:":
                        flag = True
                    if line.strip() == "":
                        flag = False
                    if flag:
                        data_frame[ac]["summary"] += line.strip() + " "


    seen_row = {}
    row = ["uniprotkb_canonical_ac","p_refseq_ac_best_match","refseq_protein_name", "refseq_protein_length","refseq_protein_summary"]
    print "\"%s\"" % ("\",\"".join(row))
    for rec in SeqIO.parse(in_file, "genbank"):
        refseq_ac_list = [rec.id]
        refseq_ac_list += rec.annotations["accessions"]
        for ac in list(set(refseq_ac_list)):
            ac_orig = ac
            ac_one = ac.split(".")[0]
            if ac.find(".") != -1:
                ac = ".".join(ac.split(".")[:-1])
            if ac not in refseq2canon:
                continue
            for feat in rec.features:
                if feat.type in ["Protein"]:
                    if ac not in data_frame:
                        continue
                    summary = data_frame[ac]["summary"]
                    summary = summary.replace("\"", "`")
                    summary = summary.replace("Summary: ", "")
                    product = feat.qualifiers["product"][0]
                    product = product.replace("\"", "`")
                    seq_len = len(str(rec.seq))
                    uniprotkb_canonical_ac = refseq2canon[ac]
                    row = [uniprotkb_canonical_ac, ac, product, str(seq_len), summary]
                    row_str = json.dumps(row)
                    if row_str not in seen_row:
                        print "\"%s\""  % ("\",\"".join(row))
                        seen_row[row_str] = True

    return



def extract_genenames_refseq_ds(species):

    seen_row = {}
    refseq2canon = {}
    canon2refseq = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        refseq = row[f_list.index("xref_id")]
        if canon not in canon2refseq:
            canon2refseq[canon] = refseq
            refseq2canon[refseq] = canon

    row = ["uniprotkb_canonical_ac","refseq_ac", "refseq_gene_name"]
    print "\"%s\""  % ("\",\"".join(row))

    seen_row = {}
    tax_id = species_obj[species]["tax_id"]
    in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_%s.gpff" % (tax_id)
    if species in ["sarscov1", "sarscov2", "hcv1", "hcv2"]:
        in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_viral.gpff"
    for record in SeqIO.parse(in_file, "genbank"):
        refseq_ac_list = [record.id]
        refseq_ac_list += record.annotations["accessions"]
        for refseq in list(set(refseq_ac_list)):
            if refseq.find(".") != -1:
                refseq = ".".join(refseq.split(".")[:-1])
            if refseq not in refseq2canon:
                continue
            canon = refseq2canon[refseq]
            for f in record.features:
                if f.type == "CDS":
                    name_list = []
                    for q in ["gene", "gene_synonym"]:
                        if q in f.qualifiers:
                            for v in f.qualifiers[q]:
                                name_list.append(v)
                    for refseq_gene_name in name_list:
                        row = [canon, refseq, refseq_gene_name]
                        row_str = json.dumps(row)
                        if row_str not in seen_row:
                            seen_row[row_str] = True
                            print "\"%s\""  % ("\",\"".join(row))
    return

def extract_ncbi_linkouts_ds(species):
    

    newrow = ["ProviderId","Database","UID","URL","IconUrl","UrlName","SubjectType","Attribute"]
    print "%s"  % (",".join(newrow))

    in_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    seen = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        refseq_ac = row[f_list.index("xref_id")]
        glygen_url = "https://glygen.org/protein/%s" % (canon.split("-")[0])
        if refseq_ac not in seen:
            newrow = ["10227","Protein",refseq_ac,glygen_url,"","","",""]
            print "%s"  % (",".join(newrow))
            seen[refseq_ac] = True
    return





def load_refseq2canon(species):

    tmp_dict_one, tmp_dict_two = {}, {}
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        refseq_isoform_ac = row[f_list.index("xref_id")]
        refseq_ac = refseq_isoform_ac.split(".")[0]
        tmp_dict_one[refseq_ac] = canon
        tmp_dict_two[refseq_isoform_ac] = canon


    return tmp_dict_one, tmp_dict_two




def extract_proteinnames_refseq_ds(species):



    refseq_ac2canon, refseq_isoform2canon = load_refseq2canon(species)
    
    refseqac2namelist = {}    
    tax_id = species_obj[species]["tax_id"]
    in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_%s.gpff" % (tax_id)
    if species in ["sarscov1", "sarscov2", "hcv1", "hcv2"]:
        in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_viral.gpff"
          
    for record in SeqIO.parse(in_file, "genbank"):
        refseq_ac_list = [record.id]
        refseq_ac_list += record.annotations["accessions"]
        for refseq_ac in list(set(refseq_ac_list)):
            if refseq_ac.find(".") != -1:
                refseq_ac = ".".join(refseq_ac.split(".")[:-1])
            if refseq_ac not in refseq_ac2canon:
                continue
            if refseq_ac not in refseqac2namelist:
                refseqac2namelist[refseq_ac] = []
            for f in record.features:
                if f.type == "Protein":
                    for q in ["product", "note"]:
                        if q in f.qualifiers:
                            for v in f.qualifiers[q]:
                                refseqac2namelist[refseq_ac].append(v)
           
    row = ["uniprotkb_canonical_ac","refseq_ac", "refseq_protein_name"]
    print "\"%s\""  % ("\",\"".join(row))

    seen_row = {}
    for refseq_isoform_ac in refseq_isoform2canon:
        refseq_ac = refseq_isoform_ac.split(".")[0]
        if refseq_ac in refseqac2namelist:
            if refseqac2namelist[refseq_ac] != []:
                for refseq_name in refseqac2namelist[refseq_ac]:
                    canon = refseq_isoform2canon[refseq_isoform_ac]
                    newrow = [canon, refseq_isoform_ac, refseq_name]
                    newrow_str = json.dumps(newrow)
                    if newrow_str not in seen_row:
                        print "\"%s\""  % ("\",\"".join(newrow))
                        seen_row[newrow_str] = True

    return





def extract_function_refseq_ds(species, tax_id):


    refseq2canon = {}
    canon2refseq = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_xref_refseq.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        refseq = row[f_list.index("xref_id")]
        canon2refseq[canon] = refseq
        refseq2canon[refseq] = canon
 

    in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_%s.gpff" % (tax_id)
    if species in ["sarscov1", "sarscov2", "hcv1", "hcv2"]:
        in_file = path_obj["downloads"] + "ncbi/refseq/current/refseq_protein_all_viral.gpff"
    
    #in_file = "tmp/toy.gpff"

    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")

    data_frame = {}
    field = "xxx"
    summary_flag = False
    with open(in_file, "r") as FR:
        for line in FR:
            newfield = line[0:12].strip() 
            field = newfield if len(newfield) > 0 else field
            value = line[12:].strip()
            if field == "VERSION":
                #ac = line.split(" ")[-1].strip()
                ac = ".".join(line.split(" ")[-1].strip().split(".")[:-1])
                data_frame[ac] = {"ac":ac, "references":[], "comment_summary":""}
            elif field == "REFERENCE":
                data_frame[ac]["references"].append({"pubmed":"", "remark":""})
            elif field == "REMARK":
                data_frame[ac]["references"][-1]["remark"] += value.strip() + " "
            elif field == "PUBMED":
                data_frame[ac]["references"][-1]["pubmed"] = value
            elif field == "COMMENT":
                if value.strip().find("Summary:") == 0:
                    summary_flag = True
                    value = value.replace("Summary:", "")
                elif value.strip() == "":
                    summary_flag = False
                if summary_flag == True:
                    data_frame[ac]["comment_summary"] += value.strip() + " "


    row = ["uniprotkb_canonical_ac","refseq_ac", "xref_key", "xref_id", "annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    seen_row = {}
    for ac in data_frame:
        comment_summary = data_frame[ac]["comment_summary"].strip()
        if comment_summary != "" and ac in refseq2canon:
            uniprotkb_canonical_ac = refseq2canon[ac]
            ann = comment_summary
            row = [uniprotkb_canonical_ac, ac, "protein_xref_refseq", ac, ann]
            row_str = json.dumps(row)
            if row_str not in seen_row:
                seen_row[row_str] = True
                print "\"%s\"" % ("\",\"".join(row))
        for o in data_frame[ac]["references"]:
            if o["remark"].find("GeneRIF:") != -1:
                o["remark"] = o["remark"][8:].strip()
                o["remark"] = o["remark"][0].upper() + o["remark"][1:]
                o["remark"] = o["remark"].replace("\"", "`")
                ann, pmid = o["remark"], o["pubmed"]
                uniprotkb_canonical_ac = "N/A"
                if ac in refseq2canon:
                    uniprotkb_canonical_ac = refseq2canon[ac]
                    row = [uniprotkb_canonical_ac, ac, "protein_xref_refseq",ac,ann]
                    row_str = json.dumps(row)
                    if row_str not in seen_row:
                        seen_row[row_str] = True
                        print "\"%s\"" % ("\",\"".join(row))
                    row = [uniprotkb_canonical_ac, ac, "protein_xref_pubmed",pmid,ann]
                    row_str = json.dumps(row)
                    if row_str not in seen_row:
                        seen_row[row_str] = True
                        print "\"%s\"" % ("\",\"".join(row))

    return




def extract_pathways_reactome_ds(species):

    row_list = [] 
    sparqlutil.load_pathways_reactome(row_list, species)

    row = ["pathway_id","pathway_name", "pathway_summary", "reaction_id_list"]
    print "\"%s\""  % ("\",\"".join(row))

    for row in row_list:
        summary = ""
        for line in row[2].split("\n"):
            if line.strip() != "":
                summary += line.strip() + " "
        row[2] = summary.strip()
        print "\"%s\""  % ("\",\"".join(row))

    return


def extract_participants_rhea_ds(species):

    seen_row = {}

    tmprow = ["reaction_id","participant_id", "participant_name","role", "xref_id", "xref_type"]
    print "\"%s\""  % ("\",\"".join(tmprow))

    row_list = []
    sparqlutil.load_enzymes_rhea(row_list, species)
    for row in row_list:
        if row[0] in ac2canon:
            tmprow = [row[1], ac2canon[row[0]],"", "enzyme", row[0],"UniProtKB"]
            row_str = json.dumps(tmprow)
            if row_str not in seen_row:
                print "\"%s\""  % ("\",\"".join(tmprow))
                seen_row[row_str] = True

    return

def extract_participants_reactome_ds(species):


    seen_row = {}


    in_file = path_obj["downloads"] + "chebi/current/database_accession.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    checbi2glytoucan = {}
    for row in data_frame["data"]:
        source = row[f_list.index("SOURCE")]
        if source != "GlyTouCan":
            continue
        glytoucan_ac = row[f_list.index("ACCESSION_NUMBER")].strip()
        xref_id = row[f_list.index("COMPOUND_ID")]
        xref_key = "glycan_xref_chebi"
        if xref_id not in checbi2glytoucan:
            checbi2glytoucan[xref_id] = []
        if glytoucan_ac not in checbi2glytoucan[xref_id]:
            checbi2glytoucan[xref_id].append(glytoucan_ac)


    sheet_obj = {"participants":{}}
    sparqlutil.load_participants_reactome(sheet_obj, [species])

    tmprow = ["reaction_id","participant_id", "participant_name","role", "xref_id", "xref_type"]
    print "\"%s\""  % ("\",\"".join(tmprow))
    
    for reaction_id in sheet_obj["participants"]:
        for role in ["input", "output"]:
            for o in sheet_obj["participants"][reaction_id][role]:
                name = o["name"].encode('ascii', 'ignore').decode('ascii')
                xref_type = o["xreftype"].encode('ascii', 'ignore').decode('ascii')
                xref_id = o["xrefid"].encode('ascii', 'ignore').decode('ascii')
                xref_id = xref_id.replace("\r", "").replace("\n", "")
                tmprow = [reaction_id,o["id"], name, role, xref_id, xref_type]
                print "\"%s\""  % ("\",\"".join(tmprow))
                if o["xreftype"] == "ChEBI" and o["xrefid"] in checbi2glytoucan:
                    for glytoucan_ac in checbi2glytoucan[o["xrefid"]]:
                        tmprow = [reaction_id,o["id"], name, role, glytoucan_ac, "GlyTouCan"]
                        row_str = json.dumps(tmprow)
                        if row_str not in seen_row:
                            print "\"%s\""  % ("\",\"".join(tmprow))
                            seen_row[row_str] = True


    row_list = []
    sparqlutil.load_enzymes_reactome(row_list, species)
    for row in row_list:
        if row[0] in ac2canon:
            tmprow = [row[1], ac2canon[row[0]],"", "enzyme", row[0],"UniProtKB"]
            row_str = json.dumps(tmprow)
            if row_str not in seen_row:
                print "\"%s\""  % ("\",\"".join(tmprow))
                seen_row[row_str] = True

    return



def extract_reactions_reactome_ds(species):

    sheet_obj = {"reactions":[]}
    sparqlutil.load_reactions(sheet_obj, [species])

    tmprow = ["reaction_id","reaction_name","cellular_location", "pmid", "pathway_id", "equation", "reaction_summary"]
    print "\"%s\""  % ("\",\"".join(tmprow))

    for o in sheet_obj["reactions"]:
        for pmid in o["pmid"].split(","):
            if o["source"] in ["Rhea"]:
                continue
            tmprow = [
                o["reactionid"],
                o["reactionname"].encode('ascii', 'ignore').decode('ascii'),
                o["cellularlocation"].encode('ascii', 'ignore').decode('ascii'),
                pmid,
                o["pathwayid"],
                o["equation"],
                o["summary"].replace("\n", "; ")
            ]
            print "\"%s\""  % ("\",\"".join(tmprow))

    return





def extract_reactions_rhea_ds(species):

    sheet_obj = {"reactions":[]}
    sparqlutil.load_reactions(sheet_obj, [species])

    tmprow = ["reaction_id","reaction_name","cellular_location", "pmid", "pathway_id", "equation", "reaction_summary"]
    print "\"%s\""  % ("\",\"".join(tmprow))

    for o in sheet_obj["reactions"]:
        if o["source"] not in ["Rhea"]:
            continue
        for pmid in o["pmid"].split(","):
            tmprow = [
                o["reactionid"],
                o["reactionname"].encode('ascii', 'ignore').decode('ascii'),
                o["cellularlocation"].encode('ascii', 'ignore').decode('ascii'),
                pmid,
                o["pathwayid"],
                o["equation"],
                o["summary"]
            ]
            print "\"%s\""  % ("\",\"".join(tmprow))

    return


def extract_misc_ds(dataset, species):

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

def extract_mutation_glycosylation_effect_ds(species, ds_name):



    fasta_file = "unreviewed/%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(fasta_file)
    n_site_info = get_glycosylation_sites(species)
    is_o_site = load_o_glyco_sites(species)
    is_glycoprotein = {}
    for combo_id in n_site_info.keys() + is_o_site.keys():
        canon = combo_id.split(" ")[0]
        is_glycoprotein[canon] = True

    ds = ds_name.replace("_glycoeffect", "")

    data_frame = {}
    in_file = path_obj["unreviewed"] +  "%s_protein_%s.csv" % (species,ds)

    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    new_f_list = ["is_known_glycoprotein", "motif_start","ref_motif", "alt_motif", "effect"]
    newrow = []
    for f in f_list + new_f_list:
        if f == "glyco_annotation":
            continue
        newrow.append(f)
    print "\"%s\""  % ("\",\"".join(newrow))


    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        aa_pos = row[f_list.index("aa_pos")]
        ref_aa = row[f_list.index("ref_aa")]
        alt_aa = row[f_list.index("alt_aa")]
        e_list = get_mutation_effect_list(seq_hash[canon], int(aa_pos),ref_aa, alt_aa)
        if e_list == []:
            continue
        newrow = []
        for f in f_list:
            if f == "glyco_annotation":
                continue
            newrow.append(row[f_list.index(f)])
        newrow.append("yes" if canon in is_glycoprotein else "no")
        for o in e_list:
            extra_row = [str(o["motifstart"]),o["refmotif"], o["altmotif"],o["effect"]]
            print "\"%s\""  % ("\",\"".join(newrow + extra_row))

    return


def extract_media_file_ds(species, ds):
   
    in_file = glob.glob("downloads/covid_media/current/%s_protein_%s*" % (species, ds))[0]
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")
    
    cmd = "cp %s unreviewed/" % (in_file)
    x = commands.getoutput(cmd)

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
    global path_obj
    global ac2canon
    global ac2canon_strict
    global canon_dict
    global genename2canon
    global canon2genenames
    global canon2geneid
    global canon2recname
    global ds2bco
    global DEBUG


    DEBUG = False
    #DEBUG = True


    species = options.species
    dataset = options.dataset




    config_obj = json.loads(open("conf/config.json", "r").read())
    path_obj = config_obj["pathinfo"]

    ds2bco = json.loads(open("generated/misc/ds2bco.json", "r").read())
    ac2canon = {}
    canon2geneid = {}
    canon2genenames, genename2canon = {}, {}
    canon2recname = {}
    canon_dict = {}
    if dataset != "masterlist":
        ac2canon, ac2canon_strict = load_ac2canon(species)
        canon2genenames, genename2canon = load_canon2genenames(species)
        canon2recname = load_canon2recname(species)
        canon_dict = load_canon_dict(species)
        canon2geneid = load_canon2geneid(species)
    
    species_obj = {}
    in_file = config_obj["pathinfo"]["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)

   
    #start logging file usage
    libgly.log_file_usage("", dataset, "write")




    data_grid = {}
    if dataset == "masterlist":
        extract_masterlist_ds(species)
    elif dataset == "transcriptlocus":
        extract_transcriptlocus_ds(species)
    elif dataset in ["recnames", "altnames", "submittednames"]:
        extract_proteinnames_ds(species, dataset)
    elif dataset == "proteinnames_refseq":
        extract_proteinnames_refseq_ds(species)
    elif dataset == "genenames_uniprotkb":
        extract_genenames_uniprotkb_ds(species)
    elif dataset == "genenames_refseq":
        extract_genenames_refseq_ds(species)
    elif dataset == "ac2pdb":
        extract_ac2pdb_ds(species)
    elif dataset == "info_uniprotkb":
        extract_info_uniprotkb_ds(species)
    elif dataset == "allsequences":
        extract_allsequences_ds(species)
    elif dataset == "canonicalsequences":
        extract_canonicalsequences_ds(species)
    elif dataset == "sequenceinfo":
        extract_sequenceinfo_ds(species)
    elif dataset == "citations_uniprotkb":
        extract_citations_uniprotkb_ds(species)
    elif dataset == "citations_reactome":
        extract_citations_reactome_ds(species)
    elif dataset in ["citations_mutation_germline", "citations_mutation_somatic",
            "citations_mutation_literature", "citations_biomarkers",
            "citations_matrixdb", "citations_iedb"]:
        extract_citations_generic_ds(species, dataset)
    elif dataset == "disease_uniprotkb":
        extract_disease_uniprotkb_ds(species)
    elif dataset == "disease_genomics_england":
        extract_disease_genomics_england_ds(species)
    elif dataset == "disease_glycosmos":
        extract_disease_glycosmos_ds(species)
    elif dataset == "disease_alliance_genome":
        extract_disease_alliance_genome_ds(species)
    elif dataset in["disease_cdg", "disease_glyco", "disease_glygen"]:
        extract_disease_compiled_ds(dataset, species)
    elif dataset == "function_uniprotkb":
        extract_function_uniprotkb_ds(species)
    elif dataset in config_obj["xref"]:
        extract_xrefs_ds(species, dataset)
    elif dataset == "go_annotation":
        extract_go_annotation_ds(species)
    elif dataset == "pro_annotation":
        extract_pro_annotation_ds(species)
    elif dataset == "site_annotation_uniprotkb":
        extract_site_annotation_uniprotkb_ds(species)
    elif dataset == "ptm_annotation_uniprotkb":
        extract_ptm_annotation_uniprotkb_ds(species)
    elif dataset == "pdb_shortlist":
        extract_pdb_shortlist_ds(species)
    elif dataset == "mutation_cancer":
        extract_mutation_cancer_ds(species)
    elif dataset == "mutation_somatic":
        extract_mutation_somatic_ds(species)
    elif dataset == "mutation_germline":
        extract_mutation_germline_ds(species)
    elif dataset == "mutation_literature":
        extract_mutation_literature_ds(species)
    elif dataset in ["mutation_cancer_glycoeffect", "mutation_somatic_glycoeffect",
            "mutation_germline_glycoeffect"]:
        extract_mutation_glycosylation_effect_ds(species, dataset)
    elif dataset == "expression_disease":
        extract_expression_disease_ds(species)
    elif dataset == "expression_normal":
        extract_expression_normal_ds(species)
    elif dataset in ["glycohydrolase", "glycosyltransferase"]:
        extract_glycoenzymes_ds(species, dataset)
    elif dataset == "genelocus":
        extract_genelocus_ds(species)
    elif dataset == "isoform_alignments":
        extract_isoform_alignments_ds(species)
    elif dataset == "info_refseq":
        extract_info_refseq_ds(species, species_obj[species]["tax_id"])
    elif dataset == "function_refseq":
        extract_function_refseq_ds(species, species_obj[species]["tax_id"])
    elif dataset == "citations_refseq":
        extract_citations_refseq_ds(species)
    elif dataset == "reactions_reactome":
        extract_reactions_reactome_ds(species)
    elif dataset == "reactions_rhea":
        extract_reactions_rhea_ds(species)
    elif dataset == "participants_reactome":
        extract_participants_reactome_ds(species)
    elif dataset == "participants_rhea":
        extract_participants_rhea_ds(species)
    elif dataset == "pathways_reactome":
        extract_pathways_reactome_ds(species)
    elif dataset == "glycosylation_motifs":
        extract_glycosylation_motifs_ds(species)
    elif dataset in ["proteinwoutsignalp", 
            "cleaved_signal_peptide_sequence", "signal_peptide_sequence"]:
        extract_misc_ds(dataset, species)
    elif dataset == "enzyme_annotation_uniprotkb":
        extract_enzyme_annotation_uniprotkb_ds(species)
    elif dataset == "binary_interactions":
        extract_binary_interactions_ds(species)
    elif dataset == "ntdata":
        extract_ntdata_ds(species)
    elif dataset == "matrixdb":
        extract_matrixdb_ds(species)
    elif dataset == "signalp_annotation":
        extract_signalp_annotation_ds(species)
    elif dataset == "integrated_site_annotation":
        extract_integrated_site_annotation_ds(species)
    elif dataset == "signalp_peptidesequences":
        extract_signalp_sequences_ds(species, "peptide_seq")
    elif dataset == "signalp_cleavedsequences":
        extract_signalp_sequences_ds(species, "cleaved_seq")
    elif dataset == "signalp_fullsequences":
        extract_signalp_sequences_ds(species, "full_seq")
    elif dataset == "ncbi_linkouts":
        extract_ncbi_linkouts_ds(species)
    elif dataset == "mutation_cancer_glycosylation_loss":
        extract_mutation_cancer_glycosylation_loss_ds(species)
    elif dataset == "biomarkers":
        extract_biomarkers_ds(species)
    elif dataset == "glycogenes":
        extract_glycogenes_ds(species)
    elif dataset == "glycogenes_glycoenzonto":
        extract_glycogenes_glycoenzonto_ds(species)
    elif dataset in ["swiss_gn2g2m3_speed2x_optimized", "swiss_gn2g2m3_speed2x"]:
        extract_media_file_ds(species, dataset)
    elif dataset == "pdb_map":
        extract_pdb_map_ds(species)
    elif dataset == "iedb":
        extract_iedb_ds(species)



    mol = "protein"
    pid = os.getpid()
    src_file = "usage/file_usage.%s.log" % (pid)
    dst_file = "usage/%s_%s_%s.fu.log" % (species, mol, dataset)

    cmd = "cat %s |sort -u > %s" % (src_file, dst_file)
    x = commands.getoutput(cmd)

    cmd = "rm -f %s" % (src_file)
    x = commands.getoutput(cmd)
    


















if __name__ == '__main__':
    
    main()


