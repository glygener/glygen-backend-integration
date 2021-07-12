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


sys.path.append('../../glytools/')
import libgly


__version__="1.0"
__status__ = "Dev"



def add_sequon_info(species, in_df):

    in_file = path_obj["unreviewed"] + "%s_protein_allsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(in_file)

    n_sequon_type_dict = {}
    in_file = "generated/misc/n_sequon_info.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    for row in data_frame["data"]:
        if row[0][-1] != "X":
            n_sequon_type_dict[row[0][-1]] = row[1]

    f_list = in_df["fields"]
    in_df["fields"] += ["n_sequon", "n_sequon_type"]
    for row in in_df["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
        g_type = row[f_list.index("glycosylation_type")]
        n_sequon, n_sequon_type = "", ""
        if aa_pos.isdigit() == True:
            aa_pos = int(aa_pos)
            if g_type.lower().find("n-linked") != -1 and canon in seq_hash and aa_pos > 0:
                n_sequon = seq_hash[canon][aa_pos-1:aa_pos+2]
                k = n_sequon[-1] 
                n_sequon_type = n_sequon_type_dict[k] if k in n_sequon_type_dict else "other"
        row += [n_sequon, n_sequon_type]

    return




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



def extract_glycosylation_unique_sources_ds(species):
   
    source_dict = {}
    file_list = glob.glob("unreviewed/%s_proteoform_glycosylation_sites_*.csv" % (species))
    for in_file in file_list:
        source = in_file.split("_sites_")[-1].split(".")[0]
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
            amino_acid = row[f_list.index("amino_acid")]
            glytoucan_ac = row[f_list.index("saccharide")]
            gly_type = row[f_list.index("glycosylation_type")].lower()
            xref_id = row[f_list.index("xref_id")]
            xref_key = row[f_list.index("xref_key")]
            src_xref_id = row[f_list.index("src_xref_id")]
            src_xref_key = row[f_list.index("src_xref_key")]

            combo_id_one = ""
            if glytoucan_ac != "":
                combo_id_one = "|".join(["1",canon, aa_pos,glytoucan_ac,gly_type,xref_key,xref_id,src_xref_key,src_xref_id])
            combo_id_two = "|".join(["2",canon, aa_pos,"",gly_type,xref_key,xref_id,src_xref_key,src_xref_id])
            combo_id_three = ""
            if xref_key == "protein_xref_pubmed":
                combo_id_three = "|".join(["3",canon,aa_pos,"",gly_type,xref_key,xref_id,src_xref_key,src_xref_id])
         

            combo_id = "1|"  + canon
            if combo_id not in source_dict:
                source_dict[combo_id] = {}
            for combo_id in [combo_id_one, combo_id_two, combo_id_three]:
                if combo_id == "":
                    continue
                if combo_id not in source_dict:
                    source_dict[combo_id] = {}
                source_dict["1|"+canon][source] = True
                source_dict[combo_id][source] = True

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","saccharide","glycosylation_type","xref_key","xref_id","src_xref_key","src_xref_id","uniqueness_flag"]
    #print "\"%s\"" % ("\",\"".join(newrow))
    out_df["fields"] = newrow

    for combo_id in sorted(source_dict):
        newrow = combo_id.split("|")
        prefix = newrow[0]
        canon = newrow[1]
        flag_name = ""
        flag_name = "unique-site-with-glycan" if prefix == "1" else flag_name
        flag_name = "unique-site" if prefix == "2" else flag_name
        flag_name = "unique-site-evidence" if prefix == "3" else flag_name
        source_list = source_dict[combo_id].keys()
        if len(source_list) == 1 and len(newrow) > 2:
            unique_src = source_list[0]
            flag = "%s-to-%s" % (flag_name, unique_src)
            if canon in source_dict:
                if source_dict[canon].keys() == [unique_src]:
                    flag += ", unique-glycoprotein-to-%s" % (source_list[0])
            newrow += [flag]
            #print "\"%s\"" % ("\",\"".join(newrow[1:]))
            out_df["data"].append(newrow[1:])
    
    return out_df


def extract_phosphorylation_citations_ds(species, ds_name):

    ds_src = ds_name.split("citations_phosphorylation_")[-1]
    black_list = get_blacklisted_pmids(species)
    data_frame = {}

    in_file = path_obj["unreviewed"] + "%s_proteoform_phosphorylation_sites_%s.csv" % (species,ds_src)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac","title","journal_name","publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    out_df["fields"] = newrow
    seen = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        src_xref_key = row[f_list.index("src_xref_key")]
        src_xref_id = row[f_list.index("src_xref_id")]
        if xref_key not in citation_xref_keys:
            continue
        if xref_id in black_list:
            continue
        newrow = libgly.get_citation(xref_id, path_obj["downloads"] + "ncbi/medline/")
        if newrow != []:
            combo_id = "%s %s" % (canon, xref_id)
            if combo_id not in seen:
                out_row = [canon] + newrow + [xref_key, xref_id,src_xref_key, src_xref_id]
                out_df["data"].append(out_row)
                seen[combo_id] = True


    return out_df

def extract_glycation_citations_ds(species, ds_name):

    ds_src = ds_name.split("citations_glycation_")[-1]
    black_list = get_blacklisted_pmids(species)
    data_frame = {}

    in_file = path_obj["unreviewed"] + "%s_proteoform_glycation_sites_%s.csv" % (species,ds_src)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac","title","journal_name","publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    out_df["fields"] = newrow
    seen = {}
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        xref_key = row[f_list.index("xref_key")]
        xref_id = row[f_list.index("xref_id")]
        src_xref_key = row[f_list.index("src_xref_key")]
        src_xref_id = row[f_list.index("src_xref_id")]
        if xref_key not in citation_xref_keys:
            continue
        if xref_id in black_list:
            continue
        newrow = libgly.get_citation(xref_id, path_obj["downloads"] + "ncbi/medline/")
        if newrow != []:
            combo_id = "%s %s" % (canon, xref_id)
            if combo_id not in seen:
                out_row = [canon] + newrow + [xref_key, xref_id,src_xref_key, src_xref_id]
                out_df["data"].append(out_row)
                seen[combo_id] = True


    return out_df


def extract_glycosylation_citations_ds(species, ds_name):

    ds_src = ds_name.split("citations_glycosylation_")[-1]
    black_list = get_blacklisted_pmids(species)
    compiled_in_file = "compiled/sarscov2_proteoform_citations_unicarbkb.csv"

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac","glytoucan_ac", "title","journal_name","publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    out_df["fields"] = newrow
    seen = {}
    log_file = "logs/%s_proteoform_%s.log" % (species, ds_name)
    FL = open(log_file, "w")
    
    in_file = path_obj["unreviewed"] + "%s_proteoform_glycosylation_sites_%s.csv" % (species,ds_src)
    
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        glytoucan_ac = row[f_list.index("saccharide")]
        xref_key = row[f_list.index("xref_key")] 
        xref_id = row[f_list.index("xref_id")]
        src_xref_key = row[f_list.index("src_xref_key")]
        src_xref_id = row[f_list.index("src_xref_id")]
        if xref_key not in citation_xref_keys:
            continue
        if xref_id in black_list:
            continue

        combo_id = "%s %s %s" % (canon, glytoucan_ac, xref_id)
        newrow = libgly.get_citation(xref_id, path_obj["downloads"] + "ncbi/medline/")
        if xref_key == "protein_xref_doi":
            newrow = libgly.get_doi_citation(xref_id, compiled_in_file)
        if newrow != []:
            if combo_id not in seen:
                out_row = [canon, glytoucan_ac] + newrow + [xref_key, xref_id, src_xref_key, src_xref_id]
                out_df["data"].append(out_row)
                #print "\"%s\"" % ("\",\"".join(out_row))
            seen[combo_id] = True
        elif xref_id not in seen:
            FL.write("%s\n" % (xref_id))
            seen[xref_id] = True
        
    FL.close()


    headers = ["uniprotkb_canonical_ac","glytoucan_ac", "title","journal_name","publication_date", "authors"]
    headers += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    in_file = "compiled/%s_proteoform_citations_unicarbkb.csv" % (species)
    if os.path.isfile(in_file) == True:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            newrow = []
            for f in headers:
                v = row[f_list.index(f)] if f in f_list else ""
                newrow.append(v)
            out_df["data"].append(newrow)
            #print "\"%s\"" % ("\",\"".join(newrow))



    return out_df

    

def load_fasta_sequences(fasta_file):
    seq_hash = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())
    return seq_hash


def load_ac2canon(species):

    is_canon = {}
    ac2canon = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        is_canon[canon] = True
        ac2canon[canon] = canon
        isoform_list = [row[f_list.index("reviewed_isoforms")], 
                        row[f_list.index("unreviewed_isoforms")]]
        for isoform in isoform_list:
            ac = isoform.split("-")[0]
            ac2canon[ac] = canon
            ac2canon[isoform] = canon

    return ac2canon, is_canon


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


def load_aa_format_dict():
    
    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_format_dict = {"one":{}, "three":{}, "glytype":{}}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        gly_type = row[f_list.index("gly_type")]
        
        form_list = [one, one.lower(), one.upper(), three, three.lower(), three.upper()]
        for f in form_list:
            aa_format_dict["three"][f] = three
            aa_format_dict["one"][f] = one
            aa_format_dict["glytype"][f] = gly_type

    return aa_format_dict


def extract_glycosylation_sites_gptwiki_ds(species):


    out_df = {"fields":[], "data":[]}
    newrow = required_fields + ["glycopeptide_id", "composition"]
    #print "\"%s\"" % ("\",\"".join(newrow))
    out_df["fields"] = newrow

    FL1 = open("logs/%s_proteoform_glycosylation_sites_gptwiki.1.log" % (species), "w")
    FL1.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    FL2 = open("logs/%s_proteoform_glycosylation_sites_gptwiki.2.log" % (species), "w")
    glycan_list = load_glycan_list()


    logged_glycan = {}
    data_frame = {}
    in_file = "downloads/gptwiki/current/glycosites.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        uniprotkb_ac = row[f_list.index("UniProt")]
        aa_pos = row[f_list.index("Site")]
        aa_three = row[f_list.index("AminoAcid")]
        saccharide = row[f_list.index("GlyTouCan")]
        composition = row[f_list.index("Composition")]
        #glycosylation_type = row[f_list.index("GlycoType")]
        glycosylation_type = "n-linked"
        glycopeptide_idlist = row[f_list.index("Glycopeptides")].split("|")
        if uniprotkb_ac not in ac2canon:
            continue
        canon = ac2canon[uniprotkb_ac]
        for glycopeptide_id in glycopeptide_idlist:
            pair_list = [["protein_xref_gptwiki", uniprotkb_ac]]
            for pair in pair_list:
                #If saccharide is not in glycan_list, force it to be ""
                if saccharide not in logged_glycan and saccharide.strip() != "" and saccharide not in glycan_list:
                    logged_glycan[saccharide] = True
                    FL2.write("\"%s\"\n" % (saccharide))
                    saccharide = ""
                g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                        pair[0],pair[1],pair[0],pair[1],glycopeptide_id,composition
                ]
                flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)
                if flag_list != []:
                    FL1.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                    continue
                out_df["data"].append(newrow)
                #print "\"%s\"" % ("\",\"".join(newrow))
    FL1.close()
    FL2.close()
    
    
    return out_df



def extract_glycosylation_sites_glyconnect_ds(species):


    in_file = path_obj["downloads"] + "glyconnect/current/glyconnect_%s.json" % (species)
    obj_list = json.loads(open(in_file, "r").read())["results"]
   
    extra_field_list = [
        'protein.id',
        'taxonomy.taxonomy_id', 'taxonomy.species',
        'structure.id', 'structure.glytoucan_id', 'structure.glycan_core', 
        'structure.glycan_type',
        'composition.format_numeric', 'composition.format_condensed', 
        'composition.format_byonic', 'composition.mass_monoisotopic', 
        'composition.mass', 'composition.format_glyconnect', 'composition.glytoucan_id',
        'source.tissue.id', 'source.tissue.name','source.tissue.uberon_id',
        'source.cell_line.id', 'source.cell_line.name', 'source.cell_line.cellosaurus_id',
        'source.cell_component.id','source.cell_component.go_id','source.cell_component.name'
    ]

    glycan_list = load_glycan_list()


    out_df = {"fields":[], "data":[]}
    newrow = required_fields
    for f in extra_field_list:
        newrow.append(f.replace(".", "_"))
    out_df["fields"] = newrow
    #print "\"%s\"" % ("\",\"".join(newrow))
    

    FL1 = open("logs/%s_proteoform_glycosylation_sites_glyconnect.1.log" % (species), "w")
    FL1.write("\"%s\"\n" % ("\",\"".join(newrow + ["filter_flags"])))
   
    FL2 = open("logs/%s_proteoform_glycosylation_sites_glyconnect.2.log" % (species), "w") 

    key_dict = {
        "taxonomy":{
            "taxonomy_id":{}
            ,"species":{}
        }
        ,"protein":{
            "id":{},
            "uniprots":{
                "uniprot_acc":{}
            }
        }
        ,"site": {
            "glyco_site": {}, 
            "location": {}
        }
        ,"structure":{
            "id":{},"glytoucan_id":{},"glycan_core":{},"glycan_type":{}
        }
        ,"composition":{
            "format_glyconnect":{},"format_condensed":{},"format_byonic":{},
            "format_numeric":{},"mass":{},"mass_monoisotopic":{},"glytoucan_id":{}
        }
        ,"source":{
            "tissue": {"id":{}, "name":{}, "uberon_id":{}}
            ,"cell_line": {"id":{}, "name":{}, "cellosaurus_id":{}}
            ,"cell_component":{"id":{}, "go_id":{}, "name":{}}
        }
    }

    logged_glycan = {}
    seen_row_str = {}
    for obj in obj_list:
        #print json.dumps(obj, indent=4)
        #print json.dumps(obj["source"], indent=4)
        #continue
        #sys.exit()

        value_dict = {}
        for k_1 in ["protein", "taxonomy", "site", "structure", "composition","source"]:
            if k_1 in obj:
                for k_2 in key_dict[k_1]:
                    combo_id = "%s.%s" % (k_1, k_2)
                    value_dict[combo_id] = ""
                    if k_2 in obj[k_1]:
                        o_2 = obj[k_1][k_2]
                        if type(o_2) is not dict:
                            value_dict[combo_id] = o_2
                        else:
                            for k_3 in o_2:
                                combo_id = "%s.%s.%s" % (k_1, k_2,k_3)
                                value_dict[combo_id] = o_2[k_3] if k_3 in o_2 else ""
        
        
        
        if "taxonomy.taxonomy_id" not in value_dict:
            continue
        if value_dict["taxonomy.taxonomy_id"] != str(species_obj[species]["tax_id"]):
            continue
        if "references" in obj:
            pmid_list, doi_list = [], []
            for o in obj["references"]:
                if "pmid" in o:
                    if o["pmid"] not in pmid_list:
                        pmid_list.append(str(o["pmid"]))
                if "doi" in o:
                    if o["doi"] not in doi_list:
                        doi_list.append(str(o["doi"]))
            value_dict["pmidlist"] = pmid_list
            value_dict["doilist"] = doi_list
        

        aa_pos, aa_three = "", ""
        if "site.glyco_site" in value_dict:
            aa_three = value_dict["site.glyco_site"].split("-")[0]
            aa_three = aa_three.split("_")[0]
            aa_pos = str(value_dict["site.location"])

        saccharide_list = []
        if "structure.glytoucan_id" in value_dict:
            saccharide_list.append(value_dict["structure.glytoucan_id"])
        if "composition.glytoucan_id" in value_dict:
            saccharide_list.append(value_dict["composition.glytoucan_id"])
        glycosylation_type = ""
        if "structure.glycan_type" in value_dict:
            glycosylation_type = value_dict["structure.glycan_type"].lower()
       
        seen_canon = {}
        if "uniprots" not in obj["protein"]:
            continue
        for o in obj["protein"]["uniprots"]:
            if "uniprot_acc" in o:
                ac = o["uniprot_acc"]
                if ac in ac2canon:
                    seen_canon[ac2canon[ac]] = True
        
        extra_value_list = []
        for f in extra_field_list:
            v = str(value_dict[f]) if f in value_dict else ""
            extra_value_list.append(v)

        pair_list = []
        if "protein.id" in value_dict:
            pair_list = [["protein_xref_glyconnect", str(value_dict["protein.id"])]]
        for pmid in list(set(value_dict["pmidlist"])):
            pair_list.append(["protein_xref_pubmed", pmid])
        #for doi in list(set(value_dict["doilist"])):
        #    pair_list.append(["protein_xref_doi", doi])
        
        for canon in seen_canon:
            for pair in pair_list:
                for saccharide in saccharide_list:
                    saccharide = saccharide.strip()
                    if saccharide == "":
                        continue
                    #If saccharide is not in glycan_list, force it to be ""
                    if saccharide not in logged_glycan and saccharide not in glycan_list:
                        logged_glycan[saccharide] = True
                        FL2.write("\"%s\"\n" % (saccharide))
                        saccharide = ""
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                        pair[0],pair[1], "protein_xref_glyconnect", str(value_dict["protein.id"])] + extra_value_list
                    row_str = json.dumps(newrow)
                    if row_str in seen_row_str:
                        continue
                    seen_row_str[row_str] = True
                    flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three,
                            glycosylation_type)
                    if flag_list == [] or flag_list == ["bad_aa"]:
                        out_df["data"].append(newrow)
                        #print "\"%s\"" % ("\",\"".join(newrow))
                    else:
                        FL1.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
    
    
    FL1.close()
    FL2.close()


    return out_df


def extract_glycosylation_sites_literature_mining_ds(species):


    known_site_dict = {}
    file_list = glob.glob("reviewed/*_proteoform_glycosylation_sites_*.csv")
    for in_file in file_list:
        if in_file.split(".")[-2] == "stat":
            continue
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
            xref_key = row[f_list.index("xref_key")]  
            if aa_pos != "" and xref_key in ["protein_xref_pubmed", "protein_xref_doi"]:
                combo_id = "%s|%s" % (canon, aa_pos)
                known_site_dict[combo_id] = True


    out_df = {"fields":[], "data":[]}
    newrow = required_fields
    #print "\"%s\"" % ("\",\"".join(newrow))
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_glycosylation_sites_literature_mining.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    glygen_ds_dict = {
        "human":"GLYDS000481",
        "mouse":"GLYDS000492",
        "rat":"GLYDS000493"
    }

    seen_row = {}
    row_list = []
    in_file = "downloads/lit_min/current/proteoform_glycosylation_sites_literature_mining.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        uniprotkb_ac = row[f_list.index("uniprotkb_ac")]
        aa_pos = row[f_list.index("glycosylation_site")]
        amino_acid = row[f_list.index("amino_acid")]
        aa_three = aa_format_dict["three"][amino_acid]
        evidence = row[f_list.index("evidence")]
        glycosylation_type = aa_format_dict["glytype"][aa_three].lower()
        saccharide = ""
        if uniprotkb_ac not in ac2canon:
            continue
        canon = ac2canon[uniprotkb_ac]
        aa_one = aa_format_dict["one"][aa_three]
        glygen_ds = glygen_ds_dict[species]
       

        pair_list = [["protein_xref_automatic_literature_mining", evidence],["protein_xref_pubmed", evidence]]
        for pair in pair_list:
            flag_list = []
            if glycosylation_type == "":
                flag_list.append("bad_glytype")
                g_type = ""
            else:
                g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()

            combo_id = "%s|%s" % (canon, aa_pos)
            if combo_id not in known_site_dict:
                flag_list.append("site_not_reported")

            newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                pair[0],pair[1], "protein_xref_automatic_literature_mining", evidence]
            flag_list += get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)

            if int(aa_pos) == 0 or int(aa_pos)  >= len(seq_hash[canon]):
                flag_list.append("bad_pos")
            elif aa_one  != seq_hash[canon][int(aa_pos)-1]:
                flag_list.append("aa_mismatch")        
            if flag_list != []:
                FL.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                continue
            newrowstr = json.dumps(newrow)
            if newrowstr not in seen_row:
                out_df["data"].append(newrow)
                #print "\"%s\"" % ("\",\"".join(newrow))
                seen_row[newrowstr] = True
    FL.close()

    return out_df


               

def extract_glycosylation_sites_literature_ds(species):


    glycan_list = load_glycan_list()

    extra_fields = []
    file_list = glob.glob("compiled/glycosylation_sites_lit_*")
    file_list += glob.glob("compiled/*_proteoform_glycosylation_sites_literature.csv")
    selected_file_list = [] 
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        canon_list = []
        for row in data_frame["data"]:
            ac_list = []
            if "uniprotkb_ac" in f_list:
                ac_list = row[f_list.index("uniprotkb_ac")].split(";")
            elif "uniprotkb_canonical_ac" in f_list:
                ac_list = row[f_list.index("uniprotkb_canonical_ac")].split(";")
            for ac in ac_list:
                if ac in ac2canon:
                    canon_list.append(ac2canon[ac])
        if canon_list == []:
            continue
        selected_file_list.append(in_file)
        
        for f in f_list:
            if f not in required_fields and f not in extra_fields:
                extra_fields.append(f)

    out_df = {"fields":[], "data":[]}
    headers = required_fields + extra_fields
    newrow = headers
    out_df["fields"] = newrow
    #print "\"%s\"" % ("\",\"".join(newrow))


    FL1 = open("logs/%s_proteoform_glycosylation_sites_literature.1.log" % (species), "w")
    FL1.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    FL2 = open("logs/%s_proteoform_glycosylation_sites_literature.2.log" % (species), "w") 

    ds_dict = {"human":"GLYDS000143", "hcv1a":"GLYDS000335", "sarscov1":"GLYDS000510"}
    src_xref_key, src_xref_id = "protein_xref_glygen_ds", ds_dict[species]

    seen_row = {}
    row_list = []
    for in_file in selected_file_list:
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            ac_list = []
            if "uniprotkb_ac" in f_list:
                ac_list = row[f_list.index("uniprotkb_ac")].split(";")
            elif "uniprotkb_canonical_ac" in f_list:
                ac_list = row[f_list.index("uniprotkb_canonical_ac")].split(";")
            pos_list = row[f_list.index("glycosylation_site_uniprotkb")].split(";")
            xref_key = row[f_list.index("xref_key")] if "xref_key" in f_list else "protein_xref_pubmed"
            xref_id = row[f_list.index("evidence")] if "evidence" in f_list else ""
            xref_id = row[f_list.index("xref_id")] if "xref_id" in f_list else xref_id
            
            aa_three = row[f_list.index("amino_acid")]
            glycosylation_type = row[f_list.index("glycosylation_type")]
            saccharide, cell_line, cellosaurus_id, abundance_normalized = "", "", "", ""
           
            if "saccharide" in f_list:
                saccharide = row[f_list.index("saccharide")]
            
            extra_values = []
            for f in extra_fields:
                v = row[f_list.index(f)] if f in f_list else ""
                extra_values.append(v)

            for ac in ac_list:
                if ac.strip() == "":
                    continue
                for pos in pos_list:
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [ac.split("-")[0], pos, aa_three,saccharide,g_type,
                            xref_key,xref_id,src_xref_key,src_xref_id] + extra_values
                    newrow_str = ",".join(newrow)
                    if newrow_str not in seen_row:
                        row_list.append(newrow)
                        seen_row[newrow_str] = True
                        if species in ["human", "hcv1a"]:
                            newrow = [ac.split("-")[0], pos, aa_three,saccharide,g_type,
                                src_xref_key,src_xref_id,src_xref_key,src_xref_id] + extra_values
                            row_list.append(newrow)
                        if species in ["hcv1a"]:
                            newrow = [ac.split("-")[0], pos, aa_three,saccharide,g_type,
                                "protein_xref_unicarbkb", ac, src_xref_key,src_xref_id] + extra_values
                            row_list.append(newrow)
    
    logged_glycan = {}
    seen_row = {}
    for row in row_list:
        uniprotkb_ac = row[0]
        if uniprotkb_ac not in ac2canon:
            continue
        canon = ac2canon[uniprotkb_ac]
        aa_pos, aa_three, glycosylation_type = row[1], row[2], row[4]
        newrow = [canon] + row[1:]
        flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)
        if flag_list != []: 
            FL1.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
            continue
        newrow_str = ",".join(newrow)
        if newrow_str in seen_row:
            continue
        seen_row[newrow_str] = True
        
        #If saccharide is not in glycan_list, force it to be ""
        saccharide = newrow[3].strip()
        if saccharide not in logged_glycan and saccharide not in glycan_list:
            logged_glycan[saccharide] = True
            FL2.write("\"%s\"\n" % (saccharide))
            newrow[3] = ""
        out_df["data"].append(newrow)
        #print "\"%s\"" % ("\",\"".join(newrow))

    FL1.close()
    FL2.close()


    return out_df


def extract_glycosylation_sites_o_glcnac_mcw_ds(species):



    extra_fields_two = ["eco_id", "carb_name", "glycosylation_subtype", "status","uniprotkb_id","gene_name",
            "recommended_name_full", "peptide"
    ]

    extra_dict_two = {}
    for ds in ["masterlist", "info_uniprotkb", "recnames"]:
        data_frame = {}
        in_file = "unreviewed/%s_protein_%s.csv" % (species, ds)
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("uniprotkb_canonical_ac")]
            if canon not in extra_dict_two:
                extra_dict_two[canon] = {}
            extra_dict_two[canon]["eco_id"] = "ECO_0000269"
            extra_dict_two[canon]["carb_name"] = "GlcNac"
            extra_dict_two[canon]["glycosylation_subtype"] = "O-GlcNAcylation"

            for f in f_list:
                if f in extra_fields_two and f not in extra_dict_two[canon]:
                    extra_dict_two[canon][f] = row[f_list.index(f)]

    
    glycan_list = load_glycan_list()
    data_frame = {}
    in_file = "downloads/mcw_oglcnac/current/%s_o-glcnacome_mcw.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
   
    extra_fields_one = []
    for f in f_list:
        if f.find("PMID #") != -1:
            continue
        if f.lower() in ["entry names","protein names","entry id","status","gene names"]:
            continue
        extra_fields_one.append(f)

    
    newrow = []
    for v in required_fields + extra_fields_one + extra_fields_two:
        v = v.strip().replace(" ", "_").replace("(", "_").replace(")", "_").lower()
        v = v.replace("__", "_").replace(",","").replace("-","_").replace("/","_")
        if v[-1] == "_":
            v = v[:-1]
        newrow.append(v)

    out_df = {"fields":[], "data":[]}
    out_df["fields"] = newrow
    #print "\"%s\"" % ("\",\"".join(newrow))
    FL = open("logs/%s_proteoform_glycosylation_sites_o_glcnac_mcw.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))


    row_list = {}
    seen = {}
    seen_row = {}
    for row in data_frame["data"]:
        flag_list = []
        ac = row[f_list.index("ENTRY ID")]
        if row[f_list.index("O-GlcNAc Sites")].strip() == "":
            flag_list.append("no o-site")
        if ac not in ac2canon:
            flag_list.append("not mapping to canonincal")
        canon = ac2canon[ac] if ac in ac2canon else ""
        
        extra_values_one = []
        for f in extra_fields_one:
            extra_values_one.append(row[f_list.index(f)])
        
        extra_values_two = []
        for f in extra_fields_two:
            v = ""
            if canon in extra_dict_two:
                v = extra_dict_two[canon][f] if f in extra_dict_two[canon] else ""
            extra_values_two.append(v)

        pmid_list = []
        for f in f_list:
            if f.find("PMID #") != -1:
                pmid = row[f_list.index(f)].strip()
                if pmid != "":
                    pmid_list.append(pmid)
        o_glcnac_score = row[f_list.index("O-GlcNAc SCORE")].strip()
        site_str = row[f_list.index("O-GlcNAc Sites")].strip()
        site_str = site_str.replace("(", "").replace(")", "")
        site_str = site_str.replace(" or ", " / ")
        site_list = site_str.split("/")
        aa_pos = ""
        aa_three = ""
        glycosylation_type = "o-linked"
        saccharide = "G49108TO"
        if saccharide not in glycan_list:
            flag_list.append("saccharide not in glycan list")

        for s in site_list:
            aa_one, aa_three, aa_pos = "","", ""
            if s.strip() != "":
                s = s.strip()
                aa_one, aa_pos = s[0], s[1:]
                aa_three = aa_format_dict["three"][aa_one]

            xref_key = "protein_xref_oglcnac_db"
            xref_id = canon.split("-")[0]
            for pmid in pmid_list:
                pair_list = [[xref_key, xref_id],["protein_xref_pubmed", pmid]]
                for pair in pair_list:
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [canon, str(aa_pos), aa_three, saccharide, g_type,
                            pair[0],pair[1], xref_key, xref_id] 
                    newrow += extra_values_one + extra_values_two
                    #relax site requirement
                    if "no o-site" in flag_list:
                        flag_list.remove("no o-site")

                    if flag_list  != []:
                        newrow[0] = ac
                        FL.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                    else:
                        flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)

                        #relax no-site rows
                        if aa_pos == "" and aa_three == "" and "aa_not_in_dictionary" in flag_list:
                            flag_list.remove("aa_not_in_dictionary")
                        if flag_list  != []:

                            FL.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                        else:
                            peptide = ""
                            if aa_pos != "":
                                aa_idx = int(aa_pos)-1
                                pep_len = 5
                                pep_one = seq_hash[canon][:aa_idx][-pep_len:]
                                pep = seq_hash[canon][aa_idx]
                                pep_two = seq_hash[canon][int(aa_pos):][:pep_len]
                                peptide = pep_one +"-"+ pep + "-" + pep_two
                                          
                            newrow[-1] = peptide
                            if canon not in row_list:
                                row_list[canon] = {}
                            if aa_pos not in row_list[canon]:
                                row_list[canon][aa_pos] = []
                            row_list[canon][aa_pos].append(newrow)

    FL.close()

    for canon in sorted(row_list.keys()):
        for aa_pos in sorted(row_list[canon].keys()):
            for newrow in row_list[canon][aa_pos]:
                out_df["data"].append(newrow)
                #print "\"%s\"" % ("\",\"".join(newrow))


    return out_df


def extract_glycosylation_sites_tyr_o_linked_ds(species):


    data_frame = {}
    in_file = "compiled/%s_proteoform_glycosylation_sites_tyr_o_linked.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    out_df = {"fields":[], "data":[]}
    newrow = required_fields
    out_df["fields"] = newrow
    #print "\"%s\"" % ("\",\"".join(newrow))
    FL = open("logs/%s_proteoform_glycosylation_sites_tyr_o_linked.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    seen = {}
    seen_row = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("uniprotkb_ac")]
        aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
        aa_three = row[f_list.index("amino_acid_uniprotkb")]
        glycosylation_type = row[f_list.index("glycosylation_type")].lower()
        evidence = row[f_list.index("evidence")]
        pmid = row[f_list.index("pmid")]
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac] if ac in ac2canon else ac
        saccharide = ""
        pair_list = [["protein_xref_glygen_ds", evidence],["protein_xref_pubmed", pmid]]
        for pair in pair_list:
            g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
            newrow = [canon, str(aa_pos), aa_three, saccharide, g_type,
                    pair[0],pair[1], "protein_xref_glygen_ds", evidence]
            flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)
            if flag_list != []:
                FL.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                continue
            newrowstr = json.dumps(newrow)
            if newrowstr not in seen_row:
                out_df["data"].append(newrow)
                #print "\"%s\"" % ("\",\"".join(newrow))
                seen_row[newrowstr] = True

    FL.close()

    return out_df





def extract_glycosylation_sites_harvard_ds(species):


    data_frame = {}
    in_file = path_obj["downloads"] + "harvard/%s_glycosylation_current.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
   
    out_df = {"fields":[], "data":[]}
    newrow = required_fields + ["composition"]
    out_df["fields"] = newrow
    #print "\"%s\"" % ("\",\"".join(newrow))
    FL1 = open("logs/%s_proteoform_glycosylation_sites_harvard.1.log" % (species), "w")
    FL1.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    FL2 = open("logs/%s_proteoform_glycosylation_sites_harvard.2.log" % (species), "w")
    glycan_list = load_glycan_list()


    logged_glycan = {}
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
                saccharide = "G70994MS"
                glycosylation_type = "o-linked"
                canon = ac2canon[ac]
                aa_one, aa_pos = s[0], s[1:]
                flag_list = []
                if aa_one not in aa_format_dict["three"]:
                    flag_list.append("invalid_aa")
                aa_three = aa_format_dict["three"][aa_one]
                flag_list += get_gly_filter_flags(species,canon,aa_pos,aa_three, 
                        glycosylation_type)
                pair_list = [["protein_xref_glygen_ds",evidence],["protein_xref_pubmed", pmid]]
                for pair in pair_list:
                    #If saccharide is not in glycan_list, force it to be ""
                    if saccharide not in logged_glycan and saccharide.strip() != "" and saccharide not in glycan_list:
                        logged_glycan[saccharide] = True
                        FL2.write("\"%s\"\n" % (saccharide))

                        saccharide = ""
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                            pair[0],pair[1],"protein_xref_glygen_ds",evidence, compo]
                    if flag_list != []:
                        FL1.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                        continue
                    newrowstr = json.dumps(newrow)
                    if newrowstr not in seen:
                        out_df["data"].append(newrow)
                        #print "\"%s\"" % ("\",\"".join(newrow))
                        seen[newrowstr] = True
    
    FL1.close()
    FL2.close()

    return out_df


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
    #in_file = path_obj["downloads"] + "glytoucan/current/export/classification.tsv"
    in_file = "unreviewed/glycan_classification.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    glytoucanac2glycosylationtype = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")].strip()
        gly_type = row[f_list.index("glycan_type")].strip().lower()
        gly_subtype = row[f_list.index("glycan_subtype")].strip()
        if glytoucan_ac not in glytoucanac2glycosylationtype:
            glytoucanac2glycosylationtype[glytoucan_ac] = []
        if gly_type not in glytoucanac2glycosylationtype[glytoucan_ac]:
            if gly_type == "n-linked":
                glytoucanac2glycosylationtype[glytoucan_ac].append(gly_type)
            if gly_type == "o-linked":
                glytoucanac2glycosylationtype[glytoucan_ac].append(gly_type)

    return glytoucanac2glycosylationtype


def load_glycan_list():

    glycan_list = []
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/glycan_properties.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_acc")]
        if ac not in glycan_list:
            glycan_list.append(ac)
    
    return glycan_list



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

    uckb2glytoucan = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/unicarbkb.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        uckb_id = row[f_list.index("UniCarbKBID")]
        if uckb_id not in uckb2glytoucan:
            uckb2glytoucan[uckb_id] = []
        uckb2glytoucan[uckb_id].append(glytoucan_ac)


    glycan_list = load_glycan_list()


    in_file = path_obj["unreviewed"] + "%s_protein_allsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(in_file)


    n_linked_aalist = []
    o_linked_aalist = []
    
    data_frame = {}
    in_file = path_obj["misc"] +  "aadict.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    aa_format_dict = {}
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        three = row[f_list.index("three")]
        one = row[f_list.index("one")]
        gly_type = row[f_list.index("gly_type")] 
        if gly_type != "":
            aa_format_dict[three] = one
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


    extra_fields_dict = {
        "sarscov2":["glytoucan_structure","swiss_model","tissue","species","abundance",
            "pdb_id","cellosaurus_id","cellosuaurus_cell_line_name","clo_id",
            "clo_cell_line_name"]
    }

    out_df = {"fields":[], "data":[]}
    newrow = required_fields + ["unicarbkb_id", "composition","curation_notes", "additonal_notes"] 
    if species in extra_fields_dict:
        newrow += extra_fields_dict[species]
    out_df["fields"] = newrow
    #print "\"%s\"" % ("\",\"".join(newrow))
    FL3.write("\"%s\"\n" % ("\",\"".join(newrow + ["qc_tag"])))


    #file_name = "all_glycosylation_current.csv"
    file_name = "unicarbkb_human_mouse_rat.csv"
    cellosaurus_info = {}
    if species in ["sarscov2"]:
        file_name = "unicarbkb_%s.csv" % (species)
        with open("downloads/cellosaurus/current/cellosaurus.txt", "r") as FR:
            for line in FR:
                if line[0:2] == "ID":
                    cl_name = line[3:].strip()
                if line[0:2] == "AC":
                    cl_id = line[3:].strip()
                    cellosaurus_info[cl_id] = {"cl_name":cl_name}

        #this is produced using make-misc*
        in_file = path_obj["unreviewed"] + "/cell_line_clo_to_cellosaurus_mapping.csv"
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            cl_id = row[f_list.index("cellosaurus_id")]
            clo_id = row[f_list.index("clo_id")]
            clo_name = row[f_list.index("cell_line_label")]
            if clo_name.find("obsolete") != -1:
                continue
            if cl_id in  cellosaurus_info:
                cellosaurus_info[cl_id]["clo_id"] = clo_id
                cellosaurus_info[cl_id]["clo_name"] = clo_name
            else:
                cellosaurus_info[cl_id] = {"clo_id":clo_id, "clo_name":clo_name}


    data_frame = {}
    in_file = path_obj["downloads"] + "unicarbkb/current/" + file_name
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    n, n1, n2 = 0, 0, 0

    for row in data_frame["data"]:
        uckb_id = row[f_list.index("Id")]
        uniprotkb_ac = row[f_list.index("Protein")].split("/")[-1]
        saccharide = row[f_list.index("Toucan")] if "Toucan" in f_list else ""
        evdn = row[f_list.index("Pmid")]
        position = row[f_list.index("Position")].split("^^")[0]
        amino_acid_type = row[f_list.index("TypeAminoAcid")]
        amino_acid = row[f_list.index("AminoAcid")].split("_")[-1]
        aa_three = ""
        if amino_acid != "":
            aa_three = amino_acid[0].upper() +  amino_acid[1:].lower()
        curation_notes = row[f_list.index("Notes")]
     
        aa_pos = 0
        if position != "" and position.isdigit():
            aa_pos = int(position)
        
        #if saccharide.find("|") == -1:
        #    continue
        #if position.find("|") == -1:
        #    continue

    
        pair_list = [["protein_xref_unicarbkb", uniprotkb_ac ]]
        if evdn not in ["", "0"]:
            pair_list.append(["protein_xref_pubmed", evdn])
        doi = ""
        if "doi" in f_list:
            doi = row[f_list.index("doi")]
            doi = doi.strip() if doi.strip() != "" else ""
            if doi.strip() != "":
                pair_list.append(["protein_xref_doi", doi])
        additional_notes = row[f_list.index("AdditionalNotes")] if "AdditionalNotes" in f_list else ""

        composition = ""
        if saccharide == "" and uckb_id[0:5] == "comp_":
            composition = uckb_id
        if uckb_id[0:5] == "comp_":
            uckb_id = ""
        extra_values = []
        if species in extra_fields_dict:
            for extra_f in extra_fields_dict[species]:
                if extra_f in f_list:
                    extra_values.append( row[f_list.index(extra_f)])
                if extra_f == "cellosaurus_id":
                    cl_id = row[f_list.index(extra_f)]
                    cl_name, clo_id, clo_name = "", "", ""
                    if cl_id in cellosaurus_info:
                        if "cl_name" in cellosaurus_info[cl_id]:
                            cl_name = cellosaurus_info[cl_id]["cl_name"]
                        if "clo_id" in cellosaurus_info[cl_id]:
                            clo_id = cellosaurus_info[cl_id]["clo_id"]
                        if "clo_name" in cellosaurus_info[cl_id]:
                            clo_name = cellosaurus_info[cl_id]["clo_name"]
                    extra_values.append(cl_name)
                    extra_values.append(clo_id)
                    extra_values.append(clo_name)


        flag_dict = {
                "no_canon":False, "bad_pos":False, "aa_mismatch":False, "invalid_aa":False,
                "glycan_without_glytype":False,"bad_glytype":False,"no_evdn":False, 
                "no_glycan_invalid_aa":False, "glycan_without_n_or_o_glytype": False
        }
        if aa_three == "":
            flag_dict["invalid_aa"] = True

        if uniprotkb_ac not in ac2canon:
            flag_dict["no_canon"] = True
            canon = uniprotkb_ac
        else:
            canon = ac2canon[uniprotkb_ac]
            if position.find("|") != -1:
                pos_list = position.split("|")
                aa_list = amino_acid.split("|")
                for j in xrange(0, len(pos_list)):
                    p = pos_list[j]
                    a = aa_list[j]
                    if p.isdigit() == False:
                        flag_dict["bad_pos"] = True
                    if int(p) == 0 or int(p) >= len(seq_hash[canon]):
                        flag_dict["bad_pos"] = True
                    if a not in aa_format_dict:
                        flag_dict["invalid_aa"] = True
                    if aa_format_dict[a] != seq_hash[canon][int(p)-1]:
                        flag_dict["aa_mismatch"] = True
            elif aa_three in aa_format_dict:
                if aa_pos == 0 or aa_pos >= len(seq_hash[canon]):
                    flag_dict["bad_pos"] = True
                elif aa_format_dict[aa_three] != seq_hash[canon][aa_pos-1]:
                    flag_dict["aa_mismatch"] = True
            else: #### amino acid is not n or o linked
                flag_dict["invalid_aa"] = True
        
        glycosylation_type = "xxx: should change by under one of the conditions"
        if saccharide == "":
            if position.find("|") != -1:
                pos_list = position.split("|")
                aa_list = amino_acid.split("|")
                if aa_list[0] in n_linked_aalist:
                    glycosylation_type = "n-linked"
                elif aa_list[0] in o_linked_aalist:
                    glycosylation_type = "o-linked"
                else:
                    glycosylation_type = ""
                    flag_dict["no_glycan_invalid_aa"] = True
            elif aa_three in n_linked_aalist:
                glycosylation_type = "n-linked"
            elif aa_three in o_linked_aalist:
                glycosylation_type = "o-linked"
            else:
                glycosylation_type = ""
                flag_dict["no_glycan_invalid_aa"] = True
        elif saccharide.find("|") != -1:
            saccharide_list = saccharide.split("|")
            tmp_list = []
            for sacc in saccharide_list:
                if sacc in glytoucanac2glycosylationtype:
                    if "n-linked" in glytoucanac2glycosylationtype[sacc]:
                        tmp_list.append("n-linked")
                    if "o-linked" in glytoucanac2glycosylationtype[sacc]:
                        tmp_list.append("o-linked")
            glycosylation_type = ";".join(tmp_list)
            if tmp_list == [] or glycosylation_type == "n-linked;o-linked":
                flag_dict["bad_glytype"] = True
        elif saccharide not in glytoucanac2glycosylationtype:
            glycosylation_type = ""
            flag_dict["glycan_without_glytype"] = True
        elif glytoucanac2glycosylationtype[saccharide] == []:
            glycosylation_type = ""
            #flag_dict["glycan_without_n_or_o_glytype"] = True
        elif position.find("|") == -1: #resume normal function if we have reported glycosylation information 
            tmp_list = []
            if saccharide in glytoucanac2glycosylationtype:
                if "n-linked" in glytoucanac2glycosylationtype[saccharide]:
                    tmp_list.append("n-linked")
                if "o-linked" in glytoucanac2glycosylationtype[saccharide]:
                    tmp_list.append("o-linked")
                if tmp_list == ["n-linked"] and aa_three not in n_linked_aalist:
                    flag_dict["n_glycan_aa_mismatch"] = True
                if tmp_list == ["o-linked"] and aa_three not in o_linked_aalist:
                    flag_dict["o_glycan_aa_mismatch"] = True
            glycosylation_type = ";".join(tmp_list)
            if glycosylation_type == "n-linked;o-linked":
                flag_dict["bad_glytype"] = True
       

        if evdn in ["", "0"] and species not in ["sarscov2"]:
            flag_dict["no_evdn"] = True
        
        if saccharide != "" and saccharide not in glycan_list:
            FL2.write("\"%s\"\n" % (saccharide))
   

        pair = pair_list[0]
        g_type = ""
        if glycosylation_type != "":
            g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
        newrow = [canon,position,aa_three,saccharide,g_type,
                pair[0],pair[1],"protein_xref_unicarbkb",
                uniprotkb_ac,uckb_id,composition,curation_notes,additional_notes
        ]
        n += 1
        flag_list = []
        for flag in flag_dict:
            if flag_dict[flag] == True:
                flag_list.append(flag)




        if flag_list != []:
            qc_tags = ";".join(flag_list)
            FL1.write("\"%s\"\n" % ("\",\"".join(newrow + [qc_tags])))
            if composition.strip() != "" and composition in uckb2glytoucan:
                for gtc_ac in uckb2glytoucan[composition]:
                    newrow[3] = gtc_ac
                    FL3.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + [qc_tags] )))
            else:
                FL3.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + [qc_tags])))
            n2 += 1
        else:
            n1 += 1
            #map composition to glytoucan (can be one-to-many)
            gtc_list = uckb2glytoucan[composition] if composition in uckb2glytoucan else []
            if composition.strip() != "" and composition in uckb2glytoucan:
                for gtc_ac in uckb2glytoucan[composition]:
                    for pair in pair_list:
                        g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                        newrow = [canon,position,aa_three,gtc_ac,g_type,
                            pair[0],pair[1],"protein_xref_unicarbkb",
                            uniprotkb_ac,uckb_id,composition,curation_notes,additional_notes
                        ]
                        FL3.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + ["validation_passed"])))
                        #print "\"%s\"" % ("\",\"".join(newrow + extra_values))
                        out_df["data"].append(newrow + extra_values)

            else:
                for pair in pair_list:
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [canon,position,aa_three,saccharide,g_type,
                        pair[0],pair[1], "protein_xref_unicarbkb",
                        uniprotkb_ac,uckb_id,composition,curation_notes,additional_notes
                    ]
                    FL3.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + ["validation_passed"])))
                    out_df["data"].append(newrow + extra_values)
                    #print "\"%s\"" % ("\",\"".join(newrow + extra_values ))
    FL1.close()
    FL2.close()
    FL3.close()

    return out_df


def extract_phosphorylation_sites_uniprotkb_ds(species):

    data_grid = {"phosphorylation":[], "genename":{}}
    sparqlutil.load_phosphorylation_sites(data_grid, species)
    sparqlutil.load_genename(data_grid, species)

    FL = open("logs/%s_proteoform_phosphorylation_sites_uniprotkb.log" % (species), "w")

    out_df = {"fields":[], "data":[]}
    row = ["uniprotkb_canonical_ac", "phosphorylation_site_uniprotkb", "amino_acid",
            "kinase_uniprotkb_canonical_ac","kinase_gene_name",
            "xref_key", "xref_id", "src_xref_key", "src_xref_id",
            "phosphorylated_residue","phosphorylation_annotation","eco_id"]
    #print "\"%s\""  % ("\",\"".join(row))
    out_df["fields"] = row
    seen_row = {}
    for row in data_grid["phosphorylation"]:
        if row[0] in ac2canon:
            ac = row[0]
            gene_name = data_grid["genename"][ac] if ac in data_grid["genename"] else ""
            xref_key, xref_id = "protein_xref_uniprotkb", ac
            src_xref_key, src_xref_id = "protein_xref_uniprotkb", ac
            canon = ac2canon[ac]
            flag_list = []
            aa_three = row[2]
            aa_one = aa_format_dict["one"][aa_three]
            aa_pos = row[3]

            if aa_pos.isdigit() == False:
                flag_list.append("bad-amino-acid")
            aa_pos = int(aa_pos) if aa_pos.isdigit() == True else -1
            if canon in seq_hash:
                if int(aa_pos) == 0 or int(aa_pos) > len(seq_hash[canon]):
                    flag_list.append("aa_pos_out_of_range")
                elif aa_one != seq_hash[canon][int(aa_pos)-1]:
                    flag_list.append("aa_mismatch")

            
            kinase_canon = ""
            kinase_gene_name = ""
            newrow = [ac2canon[ac], row[3], row[2], 
                kinase_canon, kinase_gene_name,
                xref_key, xref_id, src_xref_key, src_xref_id,
                row[1],row[4],row[7]
            ]
            if flag_list != []:
                newrow.append(";".join(flag_list))
                FL.write("\"%s\""  % ("\",\"".join(newrow)))    
                continue 
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                out_df["data"].append(newrow)
                seen_row[row_str] = True
            if row[-3] == "PubMed":
                xref_key, xref_id = "protein_xref_pubmed", row[-2]
                newrow = [ac2canon[ac], row[3], row[2],
                    kinase_canon, kinase_gene_name,
                    xref_key, xref_id, src_xref_key, src_xref_id,
                    row[1],row[4],row[7]
                ]
                row_str = json.dumps(newrow)
                if row_str not in seen_row:
                    out_df["data"].append(newrow)
                    seen_row[row_str] = True
    FL.close()

    return out_df


def extract_phosphorylation_sites_iptmnet_ds(species):
    
    canon2genename = load_canon2genename(species)

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac", "phosphorylation_site_uniprotkb", "amino_acid",
        "kinase_uniprotkb_canonical_ac","kinase_gene_name", 
        "xref_key", "xref_id", "src_xref_key", "src_xref_id"
    ]
    #print "\"%s\""  % ("\",\"".join(newrow))
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_phosphorylation_sites_iptmnet.log" % (species), "w")

    seen_row = {}
    data_frame = {}
    in_file = "downloads/iptmnet/current/ptm.txt"
    with open(in_file, "r") as FR:
        for line in FR:
            row = line.split("\t")
            flag_list = []
            if row[0] not in ["PHOSPHORYLATION"]:
                flag_list.append("no-phosphorylation-record")
            cond_list = [row[4].lower().find("homo sapiens") != -1]
            cond_list += [row[4].lower().find("mus musculus") != -1]
            cond_list += [row[4].lower().find("rattus norvegicus") != -1]
            if True not in cond_list:
                flag_list.append("not-species-of-interest")
            uniprotkb_ac = row[2]
            if uniprotkb_ac.find("-") != -1 and uniprotkb_ac not in is_canon:
                continue
            if uniprotkb_ac not in ac2canon:
                flag_list.append("uniprotkb_ac-not-mapping")
            canon = ac2canon[uniprotkb_ac] if uniprotkb_ac in ac2canon else ""
            aa_one = row[5][0]
            aa_three = aa_format_dict["three"][aa_one]
            aa_pos = row[5][1:]
            if aa_pos.isdigit() == False:
                flag_list.append("bad-amino-acid")
            aa_pos = int(aa_pos) if aa_pos.isdigit() == True else -1
            
            if canon in seq_hash:
                if int(aa_pos) == 0 or int(aa_pos) > len(seq_hash[canon]):
                    flag_list.append("aa_pos_out_of_range")
                elif aa_one != seq_hash[canon][int(aa_pos)-1]:
                    flag_list.append("aa_mismatch")
            kinase_canon = ac2canon[row[6]] if row[6] in ac2canon else ""
            kinase_gene_name = ""
            if kinase_canon in canon2genename:
                kinase_gene_name = canon2genename[kinase_canon]
            xref_key = "protein_xref_" + row[1].lower()
            if  xref_key == "protein_xref_unip":
                xref_key = "protein_xref_uniprotkb"
            xref_id = uniprotkb_ac
            src_xref_key, src_xref_id = "protein_xref_iptmnet", uniprotkb_ac
            newrow = [canon, str(aa_pos), aa_three, kinase_canon,kinase_gene_name, 
                xref_key,xref_id, src_xref_key, src_xref_id]
            if flag_list != []:
                newrow.append(";".join(flag_list))
                FL.write("\"%s\""  % ("\",\"".join(newrow)))
                continue

            #For original xref_key that comes with the data
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                out_df["data"].append(newrow)
                seen_row[row_str] = True

            #For xref_key == protein_xref_iptmnet
            xref_id = uniprotkb_ac
            xref_key = "protein_xref_iptmnet"
            newrow = [canon, str(aa_pos), aa_three, kinase_canon,kinase_gene_name,
                                    xref_key,xref_id, src_xref_key, src_xref_id]
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                out_df["data"].append(newrow)
                seen_row[row_str] = True


            #For xref_key == protein_xref_pubmed
            if row[9].strip() != "":
                pmid_list = row[9].strip().split(",")
                xref_key = "protein_xref_pubmed"
                for xref_id in pmid_list:
                    newrow = [canon, str(aa_pos), aa_three, kinase_canon,kinase_gene_name, 
                        xref_key,xref_id, src_xref_key, src_xref_id]
                    row_str = json.dumps(newrow)
                    if row_str not in seen_row:
                        out_df["data"].append(newrow)
                        seen_row[row_str] = True

    FL.close()


    return out_df


def get_gly_filter_flags(species,canon,aa_pos, aa_three, glycosylation_type):

    
    flag_list = []
    if aa_three in aa_format_dict["one"]:
        if glycosylation_type.lower() != aa_format_dict["glytype"][aa_three]:
            flag_list.append("bad_glytype")
    if glycosylation_type.find("/") != -1:
        flag_list.append("bad_glytype")

    if aa_three.strip() == "" or aa_pos.strip() == "":
        return []
    else:
        if aa_three not in aa_format_dict["one"]:
            flag_list.append("aa_not_in_dictionary")
        elif aa_pos.isdigit() == False:
            flag_list.append("aa_pos_non_numeric")
        elif int(aa_pos) == 0 or int(aa_pos) > len(seq_hash[canon]):
            flag_list.append("aa_pos_out_of_range")
        elif aa_format_dict["one"][aa_three] != seq_hash[canon][int(aa_pos)-1]:
            flag_list.append("aa_mismatch")


    return flag_list



def extract_glycosylation_sites_uniprotkb_ds(species):

    data_grid = {"glycosylation":[]}
    sparqlutil.load_glycosylation_sites(data_grid, species)


    glysubtype_dict = {}
    data_frame = {}
    in_file = "generated/misc/glycosylation_subtype.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        prefix = row[f_list.index("comment_prefix")]
        subtype = row[f_list.index("glycosylation_subtype")]
        glysubtype_dict[prefix.strip()] = subtype

     
    out_df = {"fields":[], "data":[]}
    newrow = required_fields 
    newrow += ["uniprotkb_glycosylation_annotation_comment","data_source","eco_id",
        "uniprotkb_ftid", "carb_name", "glycosylation_subtype"]
    #print "\"%s\""  % ("\",\"".join(newrow))
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_glycosylation_sites_uniprotkb.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))
    
    seen_row = {}
    for row in data_grid["glycosylation"]:
        uniprotkb_ac = row[0]
        if uniprotkb_ac in ac2canon:
            canon = ac2canon[uniprotkb_ac]
            glycosylation_type = row[1].lower()
            carb_name = row[2]
            aa_three = row[3]
            aa_pos = row[4]
            comment = row[5]
            data_source = row[6]
            evdn = row[7]
            eco_id = row[8]
            uniprotkb_ftid = row[9]
            comment_prefix = comment

            if comment.find(")") != -1:
                comment_prefix = ")".join(comment.split(")")[:-1]) + ")"
            comment_prefix = comment_prefix.strip()

            glycosylation_subtype = ""
            if comment_prefix in glysubtype_dict:
                glycosylation_subtype = glysubtype_dict[comment_prefix]

            #Ignore glycation cases
            if glycosylation_type.lower() == "n-linked" and aa_three != "Asn":
                continue
            saccharide = ""

            #All rows point to UniProtKB
            xref_key = "protein_xref_uniprotkb_gly"
            xref_id = uniprotkb_ac
            src_xref_key,src_xref_id = "protein_xref_uniprotkb_gly", uniprotkb_ac
            g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
            newrow = [canon, str(aa_pos),aa_three, saccharide, g_type,
                xref_key,xref_id,src_xref_key,src_xref_id,
                comment,data_source,eco_id,uniprotkb_ftid, carb_name, glycosylation_subtype]
           
            flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)
            if flag_list != []:
                FL.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                continue 


            #Default output (with xref_key=protein_xref_unoprotkb_gly)
            newrow_str = json.dumps(newrow)
            if newrow_str not in seen_row:
                #print "\"%s\""  % ("\",\"".join(newrow))
                out_df["data"].append(newrow)
                seen_row[newrow_str] = True

            if data_source.lower() == "uniprotkb":
                continue
            if data_source == "PubMed" and evdn.isdigit() == False:
                continue
            #Now report additional rows with evidences
            elif evdn != "":
                xref_id = evdn
                xref_key = "protein_xref_" + data_source.lower()
                xref_key = "protein_xref_pdb4glycosylation" if xref_key == "protein_xref_pdb" else xref_key
                g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                newrow = [
                    canon, str(aa_pos), aa_three, saccharide, g_type,
                    xref_key,xref_id,src_xref_key,src_xref_id,
                    comment,data_source,eco_id,uniprotkb_ftid, carb_name, glycosylation_subtype
                ]
                newrow_str = json.dumps(newrow)
                if newrow_str not in seen_row:
                    #print "\"%s\""  % ("\",\"".join(newrow))
                    out_df["data"].append(newrow)
                    seen_row[newrow_str] = True


    FL.close()

    
    return out_df




def extract_glycation_sites_uniprotkb_ds(species):
    
    data_grid = {"glycosylation":[]}
    sparqlutil.load_glycosylation_sites(data_grid, species)

    out_df = {"fields":[], "data":[]}
    newrow = required_fields
    newrow += ["uniprotkb_glycation_annotation_comment","data_source","eco_id",
        "carb_name"]
    for i in xrange(0, len(newrow)):
        if newrow[i].find("glycosylation") != -1:
            newrow[i] = newrow[i].replace("glycosylation", "glycation")
    #print "\"%s\""  % ("\",\"".join(newrow))
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_glycation_sites_uniprotkb.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    seen_row = {}
    for row in data_grid["glycosylation"]:
        uniprotkb_ac = row[0]
        if uniprotkb_ac in ac2canon:
            canon = ac2canon[uniprotkb_ac]
            glycation_type = row[1].lower()
            carb_name = row[2]
            aa_three = row[3]
            aa_pos = row[4]
            comment = row[5]
            data_source = row[6]
            evdn = row[7]
            eco_id = row[8]
            uniprotkb_ftid = row[9]
            comment_prefix = comment
            if comment.find(")") != -1:
                comment_prefix = ")".join(comment.split(")")[:-1]) + ")"
            comment_prefix = comment_prefix.strip()


            #Ignore glycation cases
            if glycation_type.lower() != "n-linked":
                continue

            if aa_three not in ["Ile","Lys","Val", "Leu"]:
                continue

            saccharide = ""
            #All rows point to UniProtKB
            xref_key = "protein_xref_uniprotkb_gly"
            xref_id = uniprotkb_ac
            src_xref_key,src_xref_id = "protein_xref_uniprotkb_gly", uniprotkb_ac
            newrow = [canon, str(aa_pos),aa_three, saccharide, glycation_type,
                xref_key,xref_id,src_xref_key,src_xref_id,
                comment,data_source,eco_id,carb_name]
            flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycation_type)

            if flag_list != [] and flag_list != ["bad_glytype"]:
                FL.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                continue

            #Default output (with xref_key=protein_xref_unoprotkb_gly)
            newrow_str = json.dumps(newrow)
            if newrow_str not in seen_row:
                out_df["data"].append(newrow)
                #print "\"%s\""  % ("\",\"".join(newrow))
                seen_row[newrow_str] = True

            if data_source.lower() == "uniprotkb":
                continue
            if data_source == "PubMed" and evdn.isdigit() == False:
                continue
            #Now report additional rows with evidences
            elif evdn != "":
                xref_id = evdn
                xref_key = "protein_xref_" + data_source.lower()
                xref_key = "protein_xref_pdb4glycosylation" if xref_key == "protein_xref_pdb" else xref_key
                newrow = [
                    canon, str(aa_pos), aa_three, saccharide, glycation_type,
                    xref_key,xref_id,src_xref_key,src_xref_id,
                    comment,data_source,eco_id,carb_name
                ]
                newrow_str = json.dumps(newrow)
                if newrow_str not in seen_row:
                    out_df["data"].append(newrow)
                    #print "\"%s\""  % ("\",\"".join(newrow))
                    seen_row[newrow_str] = True


    FL.close()

    return out_df






def extract_glycosylation_sites_pdb_ds(species):


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


    out_df = {"fields":[], "data":[]}
    newrow = required_fields 
    newrow += ["pdb_id","amino_acid_chain","glycosylation_site_pdb","carb_name"]
    #print "\"%s\""  % ("\",\"".join(newrow))
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_glycosylation_sites_pdb.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    for in_file in glob.glob(path_obj["downloads"] + "pdb/current/*-linked-*-details"):
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, "\t")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            aathree_in_pdb = in_file.split("/")[-1].split("-")[-2]
            glycosylation_type = in_file.split("/")[-1].split("-")[0] + "-linked"
            glycosylation_type = glycosylation_type.lower()
            pdb_id = row[f_list.index("PDB-ID")]
            pdb_chain = row[f_list.index("PDB-Chain")]
            amino_acid = row[f_list.index("Amino-Acid")]
            pdb_pos = int(row[f_list.index("PDB-Position")])
            carb_name = row[f_list.index("Carb")]
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
            aa_three = aathree_in_pdb[0].upper() + aathree_in_pdb[1:].lower()
            aa_pos = str(canon_pos)
            saccharide = "" 
            evidence = pdb_id
            g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
            newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                "protein_xref_pdb4glycosylation",evidence,
                "protein_xref_pdb4glycosylation",evidence,
                pdb_id,pdb_chain,str(pdb_pos),carb_name
            ]
            flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)
            if flag_list != []:
                FL.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                continue
            out_df["data"].append(newrow)
            #print "\"%s\""  % ("\",\"".join(newrow))

    FL.close()

    return out_df




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


    mem_before = libgly.process_memory()
    start_time = time.time()

    global config_obj
    global sparql
    global graph_uri
    global prefixes
    global data_grid
    global species_obj
    global required_fields 

    species = options.species
    dataset = options.dataset

    global path_obj
    global seq_hash
    global aa_format_dict
    global ac2canon
    global is_canon
    global citation_xref_keys


    citation_xref_keys = ["protein_xref_pubmed", "protein_xref_doi"]


    config_obj = json.loads(open("conf/config.json", "r").read())
    path_obj = config_obj["pathinfo"]

    species_obj = {}
    in_file = config_obj["pathinfo"]["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)

    in_file = path_obj["unreviewed"] + "%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(in_file)
    aa_format_dict = load_aa_format_dict()

    ac2canon, is_canon = load_ac2canon(species)
    required_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key",
        "xref_id","src_xref_key","src_xref_id"]


    
    out_df = {}
    if dataset == "glycosylation_sites_unicarbkb":
        out_df = extract_glycosylation_sites_unicarbkb_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_uniprotkb":
        out_df = extract_glycosylation_sites_uniprotkb_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_pdb":
        out_df = extract_glycosylation_sites_pdb_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_harvard":
        out_df = extract_glycosylation_sites_harvard_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_literature":
        out_df = extract_glycosylation_sites_literature_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_literature_mining":
        out_df = extract_glycosylation_sites_literature_mining_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "phosphorylation_sites_uniprotkb":
        out_df = extract_phosphorylation_sites_uniprotkb_ds(species)
    elif dataset == "phosphorylation_sites_iptmnet":
        out_df = extract_phosphorylation_sites_iptmnet_ds(species)
    elif dataset.find("citations_glycosylation_") != -1:
        out_df = extract_glycosylation_citations_ds(species, dataset)
    elif dataset.find("citations_phosphorylation_") != -1:
        out_df = extract_phosphorylation_citations_ds(species, dataset)
    elif dataset.find("citations_glycation_") != -1: 
        out_df = extract_glycation_citations_ds(species, dataset)
    elif dataset == "glycosylation_sites_tyr_o_linked":
        out_df = extract_glycosylation_sites_tyr_o_linked_ds(species)
    elif dataset == "glycosylation_sites_o_glcnac_mcw":
        out_df = extract_glycosylation_sites_o_glcnac_mcw_ds(species)
    elif dataset == "glycosylation_sites_glyconnect":
        out_df = extract_glycosylation_sites_glyconnect_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_gptwiki":
        out_df = extract_glycosylation_sites_gptwiki_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_pro":
        out_df = extract_glycosylation_sites_pro_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_unique_sources":
        out_df = extract_glycosylation_unique_sources_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycation_sites_uniprotkb":
        out_df = extract_glycation_sites_uniprotkb_ds(species)

    print "\"%s\"" % ("\",\"".join(out_df["fields"]))
    for row in out_df["data"]:
        print "\"%s\"" % ("\",\"".join(row))

    log_file = "logs/usage-%s-%s.log" % (species, dataset)
    libgly.log_usage(log_file, mem_before, start_time)



if __name__ == '__main__':
        main()

