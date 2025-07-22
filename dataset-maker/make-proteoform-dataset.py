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


import libgly


__version__="1.0"
__status__ = "Dev"



def log_non_csv_file_use(in_file):
    cmd = "readlink -f " + in_file
    x = commands.getoutput(cmd)
    libgly.log_file_usage(x, "", "append")
    return



def add_sequon_info(species, in_df):


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


def extract_glycosylation_sites_embl_ds(species):



    
    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["source_tissue_id", "source_tissue_name"]
    extra_fields += ["source_cell_line_cellosaurus_id","source_cell_line_cellosaurus_name"]
    extra_fields += ["start_pos","end_pos", "start_aa","end_aa", "site_seq"]
    extra_fields += ["source_glycan_type","source_gene_name","composition","composition_mass"]
    extra_fields += ["n_sequon","n_sequon_type"]


    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields


    in_file = "generated/misc/n_sequon_info.csv"
    nsequon_dict = {}
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        nsequon_dict[row[f_list.index("n_sequon")]] = row[f_list.index("n_sequon_type")]



    cmp2gtc = {}
    in_file = "downloads/glytoucan/current/export/names.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        if row[f_list.index("Domain")] == "Byonic":
            gtc = row[f_list.index("GlyTouCanAccession")]
            composition = row[f_list.index("Name")]
            cmp2gtc[composition] = gtc
  

    seen_row = {}
    in_file = "downloads/embl/current/glygen_upload.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    idx = 0
    for row in data_frame["data"]:
        idx += 1
        flag_list = []
        uniprotkb_ac = row[f_list.index("uniprotkb_ac")]
        pos = row[f_list.index("glycosylation_site_uniprotkb")]
        aa_three =  row[f_list.index("amino_acid")]
        glycosylation_type =  row[f_list.index("glycosylation_type")]
        site_seq =  row[f_list.index("peptide")]
        source_tissue_id =  row[f_list.index("source_tissue_id")]
        source_tissue_name =  row[f_list.index("source_tissue")]
        source_cell_line_cellosaurus_id =  row[f_list.index("source_cell_line_cellosaurus_id")]
        source_cell_line_cellosaurus_name =  row[f_list.index("source_cell_line_cellosaurus_name")]
        tax_id_from_infile = row[f_list.index("taxonomy_id")]
        tax_id = str(species_obj[species]["tax_id"])

        source_gene_name = row[f_list.index("gene_name")]
        source_glycan_type =  row[f_list.index("glycan_type")]
        composition =  row[f_list.index("composition")].split(" ")[0]
        composition_mass =  row[f_list.index("glycan_mass")]

        src_xref_key = "protein_xref_data_submission" if species == "human" else "protein_xref_glygen_ds"
        ds = "%s_proteoform_glycosylation_sites_embl" % (species)
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        xref_key, xref_id = "protein_xref_doi", "10.1101/2023.09.13.557529"

        saccharide = ""
        if composition in cmp2gtc:
            saccharide = cmp2gtc[composition]

        start_pos, end_pos = str(pos), str(pos)
        aa_one = aa_format_dict["one"][aa_three] if aa_three in aa_format_dict["one"] else ""
        start_aa, end_aa = aa_one, aa_one        

        canon = ac2canon[uniprotkb_ac] if uniprotkb_ac in ac2canon else ""
        seqn, seqn_type = "", ""
        if canon != "" and canon in seq_hash:
            p = int(pos)
            if p > 1 and p < len(seq_hash[canon]) - 2:
                seqn = seq_hash[canon][p-1:p+2]
                seqnx = seqn[0] + "X" + seqn[2]
                seqn_type = nsequon_dict[seqnx] if seqnx in nsequon_dict else ""
      
        #skip if tax_id of a non empty canon is different from tax_id_from_infile
        if tax_id != tax_id_from_infile:
            continue


        newrow = [canon, str(pos), aa_three, saccharide, glycosylation_type, xref_key,xref_id]
        newrow += [src_xref_key,src_xref_id,source_tissue_id,source_tissue_name]
        newrow += [source_cell_line_cellosaurus_id,source_cell_line_cellosaurus_name]
        newrow += [start_pos,end_pos,start_aa,end_aa,site_seq]
        newrow += [source_glycan_type,source_gene_name,composition,composition_mass]
        newrow += [seqn,seqn_type]
        out_df["data"].append(newrow)


    return out_df




def extract_glycosylation_sites_predicted_isoglyp_ds(species):


    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["start_pos", "end_pos", "start_aa", "end_aa"]

    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields
    score_cutoff = 10.0
    file_list = glob.glob("downloads/isoglyp/outputdb/output.%s*.csv" % (species))
    
    for in_file in file_list:
        data_frame = {}
        libgly.load_sheet_commented(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            canon = row[f_list.index("Sequence Name")].split("|")[0]
            if canon not in is_canon:
                continue
            score = float(row[f_list.index("Max")])
            if score < score_cutoff:
                continue
            pos = row[f_list.index("Position")]
            start_pos, end_pos = pos, pos 
            aa_one = row[f_list.index("S/T")]
            aa_three = aa_format_dict["three"][aa_one]
            start_aa, end_aa = aa_one, aa_one
            saccharide, g_type = "", "O-linked"
            xref_key, xref_id = "", ""
            src_xref_key, src_xref_id = "", ""
        
            src_xref_key = "protein_xref_glygen_ds"
            ds = "%s_proteoform_glycosylation_sites_predicted_isoglyp" % (species)
            src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        
            xref_key = "protein_xref_glygen_ds"
            xref_id = src_xref_id 

            newrow = [canon, pos, aa_three, saccharide, g_type, xref_key, xref_id, src_xref_key, src_xref_id]
            newrow += [start_pos, end_pos, start_aa, end_aa]
            out_df["data"].append(newrow)


    return out_df


def extract_glycosylation_sites_platelet_ds(species):

    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["glycosylation_subtype", "site_type", "start_pos", "end_pos", "start_aa", "end_aa"]
    extra_fields += ["composition", "peptide", "peptide_start_pos", "peptide_end_pos",  "site_seq", "notes"]

    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields
    in_file = "compiled/%s_proteoform_glycosylation_sites_platelet.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        newrow = []
        for f in out_df["fields"]:
            if f not in f_list:
                print "missing field:", f
                exit()
            newrow.append(row[f_list.index(f)])
        out_df["data"].append(newrow)


    return out_df


def extract_glycosylation_sites_tablemaker_ds(species):

    tissueid2name = load_tissueid2name_map()
   
    
    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = [
        "source_tissue_name", "source_tissue_id", "source_cell_line_id", "do_id", "strain",
        "src_file_name", "src_row_idx"
    ]
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields
    
    file_list = glob.glob("downloads/tablemaker/current/TD*")
    for in_file in file_list:
        src_file_name = in_file.split("/")[-1]
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        src_row_idx = 0
        for row in data_frame["data"]:
            src_row_idx += 1
            canon, amino_acid, pos, g_type = "", "", "", ""
            gtc = row[f_list.index("GlyTouCan ID")]
            source_cell_line_id = row[f_list.index("Cell line ID")]
            pmid = row[f_list.index("Evidence")]
            tax_id = row[f_list.index("Species")]
            source_tissue_id = row[f_list.index("Tissue")].replace("UBERON_", "UBERON:")
            source_tissue_name = tissueid2name[source_tissue_id] if source_tissue_id in tissueid2name else ""
            do_id = row[f_list.index("Disease")]
            strain = row[f_list.index("Strain")]
            xref_key,xref_id = "glycan_xref_pubmed", pmid
            src_xref_key,src_xref_id = "glycan_xref_tablemaker", gtc
            if tax_id not in species_obj:
                continue
            if species_obj[tax_id]["short_name"] != species:
                continue 
            newrow = [
                canon, pos, amino_acid,gtc,g_type,xref_key,xref_id,src_xref_key,src_xref_id,
                source_tissue_name, source_tissue_id, source_cell_line_id,do_id, strain,
                src_file_name, str(src_row_idx)
            ]
            out_df["data"].append(newrow)
    
    return out_df


def extract_glycosylation_sites_glycosmos_ds(species):

    
    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = [
        "taxonomy_id", "taxonomy_species", "n_sequon", "n_sequon_type"
    ]
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields
    in_file = "downloads/glycosmos/current/filtered_glycoproteins_glycan_list.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        #['Protein Name', 'UniProt ID', 'Organism', 'Site', 'GlyTouCan IDs', 'Source DB']
        ac = row[f_list.index("UniProt ID")]
        gtc_list = row[f_list.index("GlyTouCan IDs")].strip().split(",")
        pos = row[f_list.index("Site")].strip()
        if pos.isdigit() == False:
            continue
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        gtc_list_f = []
        for gtc in gtc_list:
            if gtc not in glycan_dict:
                continue
            gtc_list_f.append(gtc)
        amino_acid, g_type = "Asn", "N-linked"
        taxonomy_id, taxonomy_species = "10090", row[f_list.index("Organism")].strip()
        n_sequon = seq_hash[canon][int(pos)-1] + "X" + seq_hash[canon][int(pos)+1]
        n_sequon_type = nsequon_dict[n_sequon] if n_sequon in nsequon_dict else ""

        xref_key,xref_id = "protein_xref_glycosmos", ac
        src_xref_key = "protein_xref_glygen_ds"
        ds = "%s_proteoform_glycosylation_sites_glycosmos" % (species)
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        for gtc in gtc_list_f:
            newrow = [
                canon, pos, amino_acid,gtc,g_type,xref_key,xref_id,src_xref_key,src_xref_id,
                taxonomy_id, taxonomy_species, n_sequon, n_sequon_type
            ]
            out_df["data"].append(newrow)

    return out_df

    

def extract_glycosylation_sites_c_man_ds(species):


    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = [
        "uniprotkb_glycosylation_annotation_comment","data_source","eco_id","uniprotkb_ftid",
        "carb_name","notes","glycosylation_subtype","start_pos","end_pos","start_aa","end_aa","site_seq"
    ]
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields
    in_file = "compiled/%s_proteoform_glycosylation_sites_c_man.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        newrow = []
        for f in out_df["fields"]:
            if f not in f_list:
                print "missing field:", f
                exit()
            newrow.append(row[f_list.index(f)])
        out_df["data"].append(newrow)


    return out_df




def extract_glycosylation_sites_platelet_ds_old(species):

    composition2gtc = {}
    in_file = "generated/misc/platelet_olinked_mapping.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        composition = row[f_list.index("Glycans NHFAGNa")]
        gtc = row[f_list.index("glytoucan_ac")]
        composition2gtc[composition] = gtc



    in_file = "downloads/user_submission/platelet_o_linked/current/input.1.tsv"
    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["glycosylation_subtype", "site_type", "start_pos", "end_pos", "start_aa", "end_aa"]
    extra_fields += ["composition", "peptide", "peptide_start_pos", "peptide_end_pos",  "site_seq", "notes"]  

 
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields

    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    aa_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for row in data_frame["data"]:
        ac = row[f_list.index("Protein Name")].split("|")[1]
        peptide = row[f_list.index("Peptide")][2:-2]
        composition = row[f_list.index("Glycans NHFAGNa")]
        g_type = row[f_list.index("Glycosylation site Localisation Assignment")]
        offset = row[f_list.index("Starting position")]
        saccharide = composition2gtc[composition] if composition in composition2gtc else ""
        if g_type.lower() == "ambiguous":
            continue
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac]
        if offset.isdigit() == False:
            continue
        offset = int(offset)  

        xref_id = "38237698"
        xref_key = "protein_xref_pubmed"
        src_xref_key = "protein_xref_data_submission" if species == "human" else "protein_xref_glygen_ds"
        ds = "%s_proteoform_glycosylation_sites_platelet" % (species)
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        source_cell_line_cellosaurus_id = "CL:0000233"
        source_cell_line_cellosaurus_name = "platelet"

        notes_one = "Thrombin-activated platelet releasate proteins were found to be enriched for a wide range of O-glycan modifications. Some C-mannosylation glycosylation sites were also identified. Mutation of O-fucosylation sites within the EMI domain of MMRN1 affects secretion: T216A reduces secretion to at least 50%, whereas mutation of T1055A almost abolishes secretion. Fucosylation of these sites is carried out either by POFUT1 or a novel POFUT, but it is not carried out by POFUT2. Fucosylation of MMRN1 at T216, may represent a new POFUT1 O-fucosylation motif (C1-X-X-X-X-T-X) that is missing the typical C-terminal cysteine residue." 
        notes_two = "Thrombin-activated platelet releasate proteins were found to be enriched for a wide range of O-glycan modifications. Some C-mannosylation glycosylation sites were also identified."

        notes = notes_one if ac == "Q13201" else notes_two

        pos_dict = {}
        pos = offset - 1
        peptide_seq = ""
        for j in range(0, len(peptide)):
            if aa_set.find(peptide[j]) != -1:
                pos += 1
                peptide_seq += peptide[j]
                #print (pos, peptide[j], peptide)
            if peptide[j] == "[":
                aa = peptide[j-1]
                if peptide[j-1] == seq_hash[canon][pos-1] and peptide[j-1] in aa_format_dict["three"]:
                    pos_dict[pos] = peptide[j-1]
        
        for pos in sorted(pos_dict):
            aa_one = pos_dict[pos]
            if aa_one not in ["S", "T"]:
                continue
            aa_three = aa_format_dict["three"][aa_one]
            start_pos, end_pos = str(pos), str(pos)
            start_aa, end_aa, site_seq = aa_one, aa_one, aa_one 
            peptide_start_pos, peptide_end_pos = str(offset), str(offset + len(peptide_seq) - 1) 
            g_type = "O-linked"
            g_subtype = "O-fucosylation" if composition[0:3] == "Fuc" else ""
            site_type = "known"
            newrow = [canon, str(pos), aa_three, saccharide, g_type, xref_key, xref_id, src_xref_key,src_xref_id]
            newrow += [g_subtype, site_type, start_pos, end_pos, start_aa, end_aa]
            newrow += [composition, peptide_seq, peptide_start_pos, peptide_end_pos, site_seq, notes] 
            out_df["data"].append(newrow)


    in_file = "compiled/%s_proteoform_glycosylation_sites_platelet.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        newrow = []
        for f in out_df["fields"]:
            newrow.append(row[f_list.index(f)])
        out_df["data"].append(newrow)

         

    return out_df



def extract_glycosylation_sites_carbbank_ds(species):


    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["start_pos", "end_pos", "start_aa", "end_aa"]
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields

    in_file = "compiled/%s_proteoform_glycosylation_sites_carbbank.csv" % (species)
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        newrow = []
        for f in out_df["fields"]:
            if f not in f_list:
                print "missing field:", f
                exit()
            newrow.append(row[f_list.index(f)])
        out_df["data"].append(newrow)


    return out_df






def extract_glycosylation_sites_carbbank_ds_old(species):

    carbbank_species_map = json.loads(open("generated/misc/carbbank_species.json", "r").read())

    carbbank2glytoucan = load_carbbank2glytoucan()

    #"CC","AU","TI","CT","SC","BS","SB","DA","FC","SI","ST","MT","TN","VR","AM","PM","AN","NT","AG","DB","NC","PA","BA"

    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["start_pos", "end_pos", "start_aa", "end_aa"]
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields

    in_file = "downloads/carbbank/current/carbbank.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        db = row[f_list.index("DB")]
        carbbank_id = row[f_list.index("CC")].split(":")[1]
        mt = row[f_list.index("MT")]
        bs = row[f_list.index("BS")]
        g_type = ""
        g_type = "N-linked" if mt in ["N-linked glycopeptide", "N-linked glycoprotein"] else g_type
        g_type = "O-linked" if mt in ["O-linked glycopeptide", "O-linked glycoprotein"] else g_type
        sp = ""
        if bs.find("(CN)") != -1:
            sp = bs.split("(CN)")[1].strip().split(",")[0]
        if sp == "":
            continue
        if g_type == "":
            continue
        if sp not in carbbank_species_map:
            continue
        if species != carbbank_species_map[sp]:
            continue

        ac, xref_id = "", ""
        if db.find("SwissProt:") != -1:
            ac = db.split("SwissProt:")[1].split(" ")[0]
        if db.find("PMID:") != -1:
            xref_id = db.split("PMID:")[1].split(" ")[0]
        if xref_id == "":
            continue 
    
        canon = ac2canon[ac] if ac in ac2canon else ""
        pos, start_pos, end_pos = "", "", ""
        start_aa, end_aa = "", ""
        amino_acid  = ""
        xref_key = "protein_xref_pubmed"
        src_xref_key = "protein_xref_glygen_ds"
        ds = "%s_proteoform_glycosylation_sites_carbbank" % (species)
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        

        if carbbank_id in carbbank2glytoucan:
            for glytoucan_ac in carbbank2glytoucan[carbbank_id]:
                newrow = [canon, pos, amino_acid, glytoucan_ac, g_type, xref_key, xref_id, src_xref_key,src_xref_id]
                newrow += [start_pos, end_pos, start_aa, end_aa]
                out_df["data"].append(newrow)


    return out_df




def extract_glycosylation_sites_diabetes_glycomic_ds(species):


    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["abundance", "source_tissue_id", "source_tissue_name"]
    extra_fields += ["start_pos","end_pos", "start_aa","end_aa"]
    extra_fields += ["sample_id"]

    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields


    fig2gtc = {}
    in_file = "generated/misc/fig_to_gtc.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        fig2gtc[row[f_list.index("ID")]] = row[f_list.index("GlyTouCan")]

    log_file = "logs/%s_proteoform_glycosylation_sites_diabetes_glycomic.log" % (species)
    FL = open(log_file, "w")

    seen_row = {}
    in_file = "downloads/zagreb/current/HG_FinnRisk.txt"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]

    FL.write("\"%s\"\n" % ("\",\"".join(f_list + ["flag_list"])))
    idx = 0
    for row in data_frame["data"]:
        idx += 1
        ab_obj = {}
        for f in f_list:
            if f.find("GP") == 0:
                ab_obj[f] = row[f_list.index(f)]
        sample_id = row[f_list.index("Sample")]
        canon = ""
        pos, aa_one, aa_three = "", "", ""
        start_pos, end_pos = pos, pos
        start_aa, end_aa = aa_one, aa_one
        glycosylation_type = "N-linked"
        source_tissue_id = "UBERON:0001969"
        source_tissue_name = "blood plasma"
        src_xref_key = "protein_xref_data_submission" if species == "human" else "protein_xref_glygen_ds"
        ds = "%s_proteoform_glycosylation_sites_diabetes_glycomic" % (species)
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        xref_key, xref_id = "protein_xref_pubmed", "28905229"
        
        tmp_row_list = []
        n_empty_abdn = 0
        n_unmapped_glycans = 0
        for fig in ab_obj:
            abdn = ab_obj[fig].replace(",", ".")
            if abdn == "":
                n_empty_abdn += 1
            saccharide = fig2gtc[fig] if fig in fig2gtc else ""
            if saccharide not in glycan_dict:
                n_unmapped_glycans += 1
            if abdn != "" and saccharide in glycan_dict:   
                newrow = [canon,pos,aa_three,saccharide,glycosylation_type,xref_key,xref_id,src_xref_key,src_xref_id]
                newrow += [abdn, source_tissue_id, source_tissue_name, start_pos, end_pos, start_aa, end_aa, sample_id]
                tmp_row_list.append(newrow)    
        n_figs = len(ab_obj.keys())
        flag_list = []
        if n_empty_abdn == n_figs:
            flag_list.append("no-abundancy-info")
        if n_unmapped_glycans == n_figs:
            flag_list.append("glycans-not-in-masterlist")

        if flag_list != []:
            flag_list = list(set(flag_list))
            FL.write("\"%s\"\n" % ("\",\"".join(row + [";".join(flag_list)])))
        else:
            for newrow in tmp_row_list:
                s = json.dumps(newrow)
                if s not in seen_row:
                    out_df["data"].append(newrow)
                    seen_row[s] = True

    FL.close()

    return out_df



def load_nsequon_dict():

    in_file = "generated/misc/n_sequon_info.csv"
    nsequon_dict = {}
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        nsequon_dict[row[f_list.index("n_sequon")]] = row[f_list.index("n_sequon_type")]

    return nsequon_dict





def extract_glycosylation_sites_pdc_ccrcc_ds(species):


    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key", "xref_id","src_xref_key","src_xref_id"]
    extra_fields = ["composition","abundance", "source_tissue_id", "source_tissue_name"]
    extra_fields += ["start_pos","end_pos", "start_aa","end_aa", "site_seq", "site_type"]
    extra_fields += ["sequon", "sequon_type", "biospecimen_id"]
 
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += extra_fields

    canon2genename, genename2canon = load_canon2genename(species)
    f_list_one = ["Gene", "Glycosite", "Nglycan", "Stripped_Sequence", "Nglycan"]

    nglycan2gtc = {}
    in_file = "generated/misc/pdc_glytoucan_mapping.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        nglycan2gtc[row[f_list.index("Nglycan")]] = row[f_list.index("glytoucan")]


    log_file = "logs/%s_proteoform_glycosylation_sites_pdc_ccrcc.log" % (species)
    FL = open(log_file, "w")

    seen_row = {}
    in_file = "downloads/pdc/current/ccRCC_TMT_intact_glycopeptide_abundance_MD-MAD.tsv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
   
 
    FL.write("\"%s\"\n" % ("\",\"".join(f_list + ["flag_list"])))

    idx = 0
    for row in data_frame["data"]:
        flag_list = []
        idx += 1
        obj = {}
        for f in f_list_one:
            obj[f] = row[f_list.index(f)]
        gene_name = obj["Gene"]
        nglycan = obj["Nglycan"]
        ab_obj = {}
        for f in f_list:
            if f.find("CPT") == 0:
                ab_obj[f] = row[f_list.index(f)]
        canon_list = []
        if gene_name not in genename2canon:
            flag_list.append("no-genename2canon-map")
        else:
            canon_list = genename2canon[gene_name].keys()
        
        pos = obj["Glycosite"].replace("N", "")
        aa_one = obj["Glycosite"][0]
        aa_three = aa_format_dict["three"][aa_one] if aa_one in aa_format_dict["three"] else ""
        if pos.isdigit() == False:
            flag_list.append("bad-position")
            pos = -1
        else:
            pos = int(pos)
        
        saccharide = nglycan2gtc[nglycan] if nglycan in nglycan2gtc else ""
        composition = nglycan
        if saccharide not in glycan_dict:
            flag_list.append("glycan-not-in-masterlist")

        site_type = "known_site"   
        glycosylation_type = "N-linked"
        start_pos, end_pos = pos, pos
        start_aa, end_aa = aa_one, aa_one
        site_seq = obj["Stripped_Sequence"]
        source_tissue_id = "UBERON:0002113"
        source_tissue_name = "kidney"

        pair_list = [["protein_xref_pubmed", "37074911"]]
        src_xref_key , src_xref_id =  "protein_xref_pdc", "PDC000471"
        tmp_row_list = []
        for canon in canon_list:
            seqn, seqnx, seqn_type = "", "", ""
            if canon in seq_hash:
                if pos > 1 and pos < len(seq_hash[canon]) - 2:
                    seqn = seq_hash[canon][pos-1:pos+2]
            if seqn == "":
                flag_list.append("pos-out-of-range")
            else:
                seqnx = seqn[0] + "X" + seqn[2]
                seqn_type = nsequon_dict[seqnx] if seqnx in nsequon_dict else ""
            for pair in pair_list:
                xref_key, xref_id = pair[0], pair[1]
                n_empty_abdn = 0
                for biospecimen_id in ab_obj:
                    abdn = ab_obj[biospecimen_id]
                    abdn = "" if abdn == "NA" else abdn
                    if abdn == "":
                        n_empty_abdn += 1
                    newrow = [canon, str(pos), aa_three, saccharide, glycosylation_type]
                    newrow += [xref_key, xref_id, src_xref_key , src_xref_id]
                    newrow += [composition,abdn,source_tissue_id, source_tissue_name]
                    newrow += [str(start_pos), str(end_pos), start_aa, end_aa, site_seq, site_type]
                    newrow += [seqn, seqn_type, biospecimen_id]
                    if flag_list == []:
                        tmp_row_list.append(newrow)
                #if n_empty_abdn == len(ab_obj.keys()):
                #    flag_list.append("empty-abdn-info")

        if flag_list != []:
            flag_list = list(set(flag_list))
            FL.write("\"%s\"\n" % ("\",\"".join(row + [";".join(flag_list)])))
        else:
            for newrow in tmp_row_list:
                s = json.dumps(newrow)
                if s not in seen_row:
                    out_df["data"].append(newrow)
                    seen_row[s] = True
    FL.close()


    return out_df
    

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
                combo_id_one = "^|^".join(["1",canon, aa_pos,glytoucan_ac,gly_type,xref_key,xref_id,src_xref_key,src_xref_id])
            combo_id_two = "^|^".join(["2",canon, aa_pos,"",gly_type,xref_key,xref_id,src_xref_key,src_xref_id])
            combo_id_three = ""
            xref_id = xref_id.replace("GLYDS", "GLY_")
            src_xref_id = src_xref_id.replace("GLYDS", "GLY_")
            if xref_key == "protein_xref_pubmed":
                combo_id_three = "^|^".join(["3",canon,aa_pos,"",gly_type,xref_key,xref_id,src_xref_key,src_xref_id])
         

            combo_id = "1^|^"  + canon
            if combo_id not in source_dict:
                source_dict[combo_id] = {}
            for combo_id in [combo_id_one, combo_id_two, combo_id_three]:
                if combo_id == "":
                    continue
                if combo_id not in source_dict:
                    source_dict[combo_id] = {}
                source_dict["1^|^"+canon][source] = True
                source_dict[combo_id][source] = True

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb","saccharide","glycosylation_type","xref_key","xref_id","src_xref_key","src_xref_id","uniqueness_flag"]
    out_df["fields"] = newrow

    for combo_id in sorted(source_dict):
        newrow = combo_id.split("^|^")
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
            out_df["data"].append(newrow[1:])
    
    return out_df


def extract_phosphorylation_citations_ds(species, ds_name):

    ds_src = ds_name.split("citations_phosphorylation_")[-1]
    black_list = get_blacklisted_pmids(species)
    data_frame = {}

    in_file = path_obj["unreviewed"] + "%s_proteoform_phosphorylation_%s.csv" % (species,ds_src)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac","title","journal_name","publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    out_df["fields"] = newrow
    seen = {}
    log_file = path_obj["logs"] + "/%s_proteoform_%s.log"  % (species, ds_name)
    FL = open(log_file, "w")
    log_dict = {}
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
        cite_info = libgly.get_citation(xref_id)
        newrow = cite_info["row"]
        if newrow != []:
            combo_id = "%s %s" % (canon, xref_id)
            if combo_id not in seen:
                out_row = [canon] + newrow + [xref_key, xref_id,src_xref_key, src_xref_id]
                out_df["data"].append(out_row)
                seen[combo_id] = True
        
        if newrow == [] and xref_id not in log_dict:
            FL.write("%s,%s\n" % (xref_id, ";".join(cite_info["flaglist"])))
            log_dict[xref_id] = True


    FL.close()

    return out_df

def extract_glycation_citations_ds(species, ds_name):

    ds_src = ds_name.split("citations_glycation_")[-1]
    black_list = get_blacklisted_pmids(species)
    data_frame = {}

    in_file = path_obj["unreviewed"] + "%s_proteoform_glycation_%s.csv" % (species,ds_src)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac","title","journal_name","publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id"]
    out_df["fields"] = newrow
    seen = {}

    log_file = path_obj["logs"] + "/%s_proteoform_%s.log"  % (species, ds_name)
    FL = open(log_file, "w")
    log_dict = {}

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
        cite_info = libgly.get_citation(xref_id)
        newrow = cite_info["row"]
        if newrow != []:
            combo_id = "%s %s" % (canon, xref_id)
            if combo_id not in seen:
                out_row = [canon] + newrow + [xref_key, xref_id,src_xref_key, src_xref_id]
                out_df["data"].append(out_row)
                seen[combo_id] = True
        if newrow == [] and xref_id not in log_dict:
            FL.write("%s,%s\n" % (xref_id, ";".join(cite_info["flaglist"])))
            log_dict[xref_id] = True

    FL.close()

    return out_df


def extract_glycosylation_citations_ds(species, ds_name):

    ds_src = ds_name.split("citations_glycosylation_")[-1]
    black_list = get_blacklisted_pmids(species)
    compiled_in_file = "compiled/doi_citations.csv"


    out_df = {"fields":[], "data":[]}

    newrow = ["uniprotkb_canonical_ac","title","journal_name","publication_date", "authors"]
    newrow += ["xref_key", "xref_id", "src_xref_key", "src_xref_id", "glytoucan_ac"]
    out_df["fields"] = newrow
    seen = {}
     
    in_file = path_obj["unreviewed"] + "%s_proteoform_glycosylation_%s.csv" % (species,ds_src)

    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    
    log_file = path_obj["logs"] + "/%s_proteoform_%s.log"  % (species, ds_name)
    FL = open(log_file, "w")
    log_dict = {}

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
        if src_xref_id == "":
            src_xref_key, src_xref_id = xref_key, xref_id

        combo_id = "%s %s %s" % (canon, glytoucan_ac, xref_id)
        cite_info = libgly.get_citation(xref_id)
        newrow = cite_info["row"]
        if xref_key == "protein_xref_doi":
            cite_info = libgly.get_doi_citation(xref_id, compiled_in_file)
            newrow = cite_info["row"]
        if newrow != []:
            if combo_id not in seen:
                out_row = [canon] + newrow + [xref_key, xref_id, src_xref_key, src_xref_id, glytoucan_ac]
                out_df["data"].append(out_row)
            seen[combo_id] = True
        if newrow == [] and xref_id not in log_dict:
            FL.write("%s,%s\n" % (xref_id, ";".join(cite_info["flaglist"])))
            log_dict[xref_id] = True

        
    FL.close()


    headers = ["uniprotkb_canonical_ac","title","journal_name","publication_date", "authors"]
    headers += ["xref_key", "xref_id", "src_xref_key", "src_xref_id", "glytoucan_ac"]
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


    return out_df

    

def load_fasta_sequences(fasta_file):
    seq_hash = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())
    return seq_hash


def load_ac2canon_strict(species):

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

    return ac2canon

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

    genename2canon = {}
    canon2genename = {}
    data_frame = {}
    in_file = path_obj["unreviewed"] + "%s_protein_masterlist.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        canon = row[f_list.index("uniprotkb_canonical_ac")]
        gene_name = row[f_list.index("gene_name")]
        canon2genename[canon] = gene_name
        if gene_name not in genename2canon:
            genename2canon[gene_name] = {}
        genename2canon[gene_name][canon] = True
    
    return canon2genename, genename2canon

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

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()


    out_df = {"fields":[], "data":[]}
    newrow = required_output_fields + ["glycopeptide_id", "composition", "glycan_xref_key","glycan_xref_id"]
    out_df["fields"] = newrow

    FL1 = open("logs/%s_proteoform_glycosylation_sites_gptwiki.1.log" % (species), "w")
    FL1.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    FL2 = open("logs/%s_proteoform_glycosylation_sites_gptwiki.2.log" % (species), "w")


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
        site_link = row[f_list.index("SiteLink")]

        #glycosylation_type = row[f_list.index("GlycoType")]
        glycosylation_type = "n-linked"
        glycopeptide_idlist = row[f_list.index("Glycopeptides")].split("|")
        if uniprotkb_ac not in ac2canon:
            continue
        canon = ac2canon[uniprotkb_ac]
        for glycopeptide_id in glycopeptide_idlist:
            pair_list = [["protein_xref_gptwiki", site_link]]
            for pair in pair_list:
                #If saccharide is not in glycan_dict, force it to be ""
                if saccharide.strip() != "" and saccharide not in glycan_dict:
                    if saccharide not in logged_glycan:
                        logged_glycan[saccharide] = True
                        FL2.write("\"%s\"\n" % (saccharide))
                    saccharide = ""
                g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                        pair[0],pair[1],pair[0],pair[1],glycopeptide_id,composition,
                        "glycan_xref_gptwiki", saccharide
                ]
                flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)
                if flag_list != []:
                    FL1.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                    continue
                out_df["data"].append(newrow)
    FL1.close()
    FL2.close()
    
    return out_df



def extract_glycosylation_sites_unicarbkb_glycomics_study_ds(species):


    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()


    data_frame = {}
    in_file = "downloads/unicarbkb/current/known_unicarbkb_glycomics_study.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    
    newrow = required_output_fields + ["unicarbkb_id", "composition","curation_notes", "additonal_notes"]
    tmp_list = []
    for f in f_list:
        if f not in newrow:
            tmp_list.append(f)
    
    extra_fields = []
    for f in tmp_list:
        if f not in ["cell_line", "cellosaurus"]:
            extra_fields.append(f)
            f = "source_" + f if f in ["tissue_id", "tissue_name"] else f
            newrow.append(f)
   
    long2short = {}
    for sp in species_obj:
        if sp.isdigit() == False:
            long2short[species_obj[sp]["long_name"]] = sp
     
    out_df = {"fields":[], "data":[]}
    out_df["fields"] = newrow

    log_file_one = path_obj["logs"] + "%s_proteoform_glycosylation_sites_unicarbkb_glycomics_study.1.log" % (species)
    FL1 = open(log_file_one, "w")
    log_file_two = path_obj["logs"] + "%s_proteoform_glycosylation_sites_unicarbkb_glycomics_study.2.log" % (species)
    FL2 = open(log_file_two, "w")

    log_file_three = path_obj["intermediate"] + "%s_proteoform_glycosylation_sites_unicarbkb_glycomics_study.csv" % (species)
    FL = open(log_file_three, "w")

    logged_glycan = {}
    
    n, n1, n2 = 0, 0, 0
    for row in data_frame["data"]:
        uniprotkb_ac = row[f_list.index("protein")]
        position = row[f_list.index("position")]
        amino_acid = row[f_list.index("amino_acid")]
        saccharide = row[f_list.index("glytoucan_structure")]
        saccharide = row[f_list.index("toucan")] if saccharide == "" else saccharide

        composition = row[f_list.index("composition1")]
        evdn = row[f_list.index("pmid")] 
        sp = row[f_list.index("species")]
        uckb_id = row[f_list.index("Id")]
        additional_notes = row[f_list.index("additional_notes")]
        curation_notes = ""
        sp_short = long2short[sp] if sp in long2short else ""  
        if sp_short != species:
            continue
        aa_three = aa_format_dict["three"][amino_acid] if amino_acid in aa_format_dict["three"] else ""
        flag_dict = {
                "no_canon":False, "bad_pos":False, "aa_mismatch":False, "invalid_aa":False,
                "glycan_without_glytype":False,"bad_glytype":False,"no_evdn":False,
                "no_glycan_invalid_aa":False, "glycan_without_n_or_o_glytype": False
        }
        
        if uniprotkb_ac == "" or uniprotkb_ac not in ac2canon:
            flag_dict["no_canon"] = True
            canon = uniprotkb_ac
        else:
            canon = ac2canon[uniprotkb_ac]
            qc_glyco_position(position, amino_acid, canon, flag_dict) #unicarbkb_glycomics_study_ds
 
        glycosylation_type = qc_glyco_type(amino_acid,saccharide,flag_dict)

        g_type = ""
        if glycosylation_type != "":
            g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
       
        extra_values = []
        for extra_f in extra_fields:
            if f not in ["cell_line", "cellosaurus"]:
                v = row[f_list.index(extra_f)]
                if f == "tissue_name":
                    tid = row[f_list.index("tissue_id")] 
                    v = tissueid2name[tid] if tid in tissueid2name else ""
                extra_values.append(v)


        pair_list = []
        if evdn not in ["", "0"]:
            pair_list.append(["protein_xref_pubmed", evdn])

        ds = "%s_proteoform_glycosylation_sites_unicarbkb_glycomics_study" % (species)
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
            
        pair = pair_list[0]
        newrow = [canon,position,aa_three,saccharide,g_type,
            pair[0],pair[1],"protein_xref_unicarbkb_ds",src_xref_id,uckb_id,
            composition,curation_notes,additional_notes
        ]
        n += 1
        flag_list = []
        for flag in flag_dict:
            if flag_dict[flag] == True:
                flag_list.append(flag)
        #if flag_list != []:
        if flag_list != [] and flag_list != ["no_canon"]:
            qc_tags = ";".join(flag_list)
            FL1.write("\"%s\"\n" % ("\",\"".join(newrow + [qc_tags])))
            if composition.strip() != "" and composition in uckb2glytoucan:
                for gtc_ac in uckb2glytoucan[composition]:
                    newrow[3] = gtc_ac
                    FL.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + [qc_tags] )))
            else:
                FL.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + [qc_tags])))
            n2 += 1
        else:
            n1 += 1
            #map composition to glytoucan (can be one-to-many)
            if composition.strip() != "" and composition in uckb2glytoucan:
                for gtc_ac in uckb2glytoucan[composition]:
                    for pair in pair_list:
                        g_type = ""
                        if glycosylation_type != "":
                            g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                        newrow = [canon,position,aa_three,gtc_ac,g_type,
                            pair[0],pair[1],"protein_xref_unicarbkb_ds",
                            src_xref_id,uckb_id,composition,curation_notes,additional_notes
                        ]
                        FL.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + ["validation_passed"])))
                        out_df["data"].append(newrow + extra_values)
            else:
                for pair in pair_list:
                    g_type = ""
                    if glycosylation_type != "":
                        g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [canon,position,aa_three,saccharide,g_type,
                        pair[0],pair[1], "protein_xref_unicarbkb_ds",
                        src_xref_id,uckb_id,composition,curation_notes,additional_notes
                    ]
                    FL.write("\"%s\"\n" % ("\",\"".join(newrow + extra_values + ["validation_passed"])))
                    out_df["data"].append(newrow + extra_values) 

    FL1.close()
    FL2.close()
    FL.close()


    #take care of rows with unmapped tissueid and clid 
    log_file = "logs/%s_proteoform_glycosylation_sites_unicarbkb_glycomics.3.log" % (species)
    FW = open(log_file, "w")
    FW.write("\"%s\"\n" % ("\",\"".join(out_df["fields"] + ["flags"])))
    row_list = []
    for newrow in out_df["data"]:
        flag_list = []
        # log unmappable tissueid 
        if "source_tissue_id" in out_df["fields"]:
            tid = newrow[out_df["fields"].index("source_tissue_id")]
            tid = "" if tid.strip() == "0" else tid
            tid = tid.replace("UBERON_", "UBERON:")
            tname = tissueid2name[tid] if tid in tissueid2name else ""
            newrow[out_df["fields"].index("source_tissue_id")] = tid
            newrow[out_df["fields"].index("source_tissue_name")] = tname
            if tid.strip() != "" and tname == "":
                flag_list.append("unmapped_source_tissue_id(%s)" %(tid))

        # log unmappable clid
        if "source_cell_line_cellosaurus_id" in out_df["fields"]:
            clid = newrow[out_df["fields"].index("source_cell_line_cellosaurus_id")]
            clid = "" if clid.strip() == "0" else clid
            cname = clid2name[clid] if clid in clid2name else ""
            newrow[out_df["fields"].index("source_cell_line_cellosaurus_name")] = cname
            if clid.strip() != "" and cname == "":
                flag_list.append("unmapped_source_cell_line_cellosaurus_id(%s)"%(clid))

        if flag_list != []:
            FW.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
        else:
            row_list.append(newrow)
    FW.close()
    out_df["data"] = row_list


    return out_df

def extract_glycosylation_sites_glyconnect_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()

    in_file = path_obj["downloads"] + "glyconnect/current/glyconnect_%s.json" % (species)
    log_non_csv_file_use(in_file)


    obj_list = json.loads(open(in_file, "r").read())["results"]
   
    extra_field_list = [
        'protein.id',
        'taxonomy.taxonomy_id', 'taxonomy.species',
        'structure.id', 'structure.glytoucan_id', 'structure.glycan_core', 
        'structure.glycan_type',
        'composition.format_numeric', 'composition.format_condensed', 
        'composition.format_byonic', 'composition.mass_monoisotopic', 
        'composition.mass', 'composition.format_glyconnect', 'composition.glytoucan_id',
        'source.tissue.name','source.tissue.uberon_id',
        'source.cell_line.name', 'source.cell_line.cellosaurus_id',
        'source.cell_component.id','source.cell_component.go_id','source.cell_component.name'
    ]



    out_df = {"fields":[], "data":[]}
    newrow = required_output_fields
    for f in extra_field_list:
        f_new = f.replace(".", "_")
        f_new = f_new.replace("source_cell_line_name", "source_cell_line_cellosaurus_name")
        f_new = f_new.replace("source_tissue_uberon_id", "source_tissue_id")
        newrow.append(f_new)
    out_df["fields"] = newrow
    

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
        #because yeast tax_id in glyconnect is 4932, not 559292
        tmp_map = {"559292":"4932"} 
        t_id = str(species_obj[species]["tax_id"])
        t_id = tmp_map[t_id] if t_id in tmp_map else t_id
        if value_dict["taxonomy.taxonomy_id"] != t_id:
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
            gtc = value_dict["structure.glytoucan_id"]
            #If saccharide is not in glycan_dict, force it to be ""
            if gtc not in glycan_dict:
                value_dict["structure.glytoucan_id"] = ""
                saccharide_list.append("")
                if gtc not in logged_glycan:
                    logged_glycan[gtc] = True
                    FL2.write("\"%s\"\n" % (gtc))
            else:
                saccharide_list.append(gtc)
        if "composition.glytoucan_id" in value_dict:
            gtc = value_dict["composition.glytoucan_id"]
            #If saccharide is not in glycan_dict, force it to be ""
            if gtc not in glycan_dict:
                value_dict["composition.glytoucan_id"] = ""
                saccharide_list.append("")
                if gtc not in logged_glycan:
                    logged_glycan[gtc] = True
                    FL2.write("\"%s\"\n" % (gtc))
            else:
                saccharide_list.append(gtc)

        saccharide_list = list(set(saccharide_list))

        glycosylation_type = ""
        if "structure.glycan_type" in value_dict:
            glycosylation_type = value_dict["structure.glycan_type"].lower()
       
        seen_canon = {}
        if "uniprots" not in obj["protein"]:
            seen_canon[""] = True
        else:
            for o in obj["protein"]["uniprots"]:
                if "uniprot_acc" in o:
                    ac = o["uniprot_acc"]
                    if ac in ac2canon:
                        seen_canon[ac2canon[ac]] = True

        extra_value_list = []
        for f in extra_field_list:
            v = str(value_dict[f]) if f in value_dict else ""
            #if f == "source.cell_line.cellosaurus_id":
            #    v = v.replace("_", ":")
            extra_value_list.append(v)


        pair_list = []
        #if "protein.id" in value_dict:
        #    pair_list = []
        #    if str(value_dict["protein.id"]) != "None":
        #        pair_list.append(["protein_xref_glyconnect", str(value_dict["protein.id"])])
        #    elif str(value_dict["structure.id"]) != "None":
        #        pair_list.append(["glycan_xref_glyconnect", str(value_dict["structure.id"])])

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
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    src_xref_key, src_xref_id = "protein_xref_glyconnect", str(value_dict["protein.id"])
                    if src_xref_id == "None":
                        src_xref_key, src_xref_id ="glycan_xref_glyconnect",str(value_dict["structure.id"])
                    newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                        pair[0],pair[1], src_xref_key, src_xref_id] + extra_value_list
                    row_str = json.dumps(newrow)
                    if row_str in seen_row_str:
                        continue
                    seen_row_str[row_str] = True
                    flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three,
                            glycosylation_type)
                    
                    if flag_list == [] or flag_list == ["bad_aa"]:
                        out_df["data"].append(newrow)
                    else:
                        FL1.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
    
    
    FL1.close()
    FL2.close()

    #take care of rows with unmapped tissueid and clid 
    log_file = "logs/%s_proteoform_glycosylation_sites_glyconnect.3.log" % (species)
    FW = open(log_file, "w")
    FW.write("\"%s\"\n" % ("\",\"".join(out_df["fields"] + ["flags"])))
    row_list = []
    for newrow in out_df["data"]:
        flag_list = []
        # log unmappable tissueid 
        if "source_tissue_id" in out_df["fields"]:
            tid = newrow[out_df["fields"].index("source_tissue_id")]
            tid = "" if tid.strip() == "0" else tid
            tid = tid.replace("UBERON_", "UBERON:")
            tname = tissueid2name[tid] if tid in tissueid2name else ""
            newrow[out_df["fields"].index("source_tissue_id")] = tid
            newrow[out_df["fields"].index("source_tissue_name")] = tname
            if tid.strip() != "" and tname == "":
                flag_list.append("unmapped_source_tissue_id(%s)" %(tid))

        # log unmappable clid
        if "source_cell_line_cellosaurus_id" in out_df["fields"]:
            clid = newrow[out_df["fields"].index("source_cell_line_cellosaurus_id")]
            clid = "" if clid.strip() == "0" else clid
            cname = clid2name[clid] if clid in clid2name else ""
            newrow[out_df["fields"].index("source_cell_line_cellosaurus_name")] = cname
            if clid.strip() != "" and cname == "":
                flag_list.append("unmapped_source_cell_line_cellosaurus_id(%s)"%(clid))

        if flag_list != []:
            FW.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
        else:
            row_list.append(newrow)
    FW.close()
    out_df["data"] = row_list


    return out_df


def extract_glycosylation_sites_literature_mining_manually_verified_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()


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
    newrow = required_output_fields + ["mining_tool"]
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_glycosylation_sites_literature_mining_manually_verified.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    seen_row = {}
    row_list = []
    in_file = "compiled/proteoform_glycosylation_sites_literature_mining_manually_verified.csv"

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
        ds = "%s_proteoform_glycosylation_sites_literature_mining_manually_verified" % (species)
        src_xref_key = "protein_xref_glygen_ds"
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        
        pair_list = [["protein_xref_pubmed", evidence]]

        for pair in pair_list:
            flag_list = []
            if glycosylation_type == "":
                flag_list.append("bad_glytype")
                g_type = ""
            else:
                g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()

            combo_id = "%s|%s" % (canon, aa_pos)
            #if combo_id not in known_site_dict:
            #    flag_list.append("site_not_reported")

            mining_tool = "RLIMS-G"
            newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                pair[0],pair[1], src_xref_key, src_xref_id, mining_tool]
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
                seen_row[newrowstr] = True
    FL.close()

  
    return out_df


def get_known_site_dict(species):

    tmp_dict = {}
    file_list = glob.glob("reviewed/%s_proteoform_glycosylation_sites_*.csv" % (species))
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
                tmp_dict[combo_id] = True

    return tmp_dict




def extract_glycosylation_sites_literature_mining_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()


    known_site_dict = {}
    #known_site_dict = get_known_site_dict(species)

    out_df = {"fields":[], "data":[]}
    newrow = required_output_fields
    out_df["fields"] = newrow + ["mining_tools"]

    FL = open("logs/%s_proteoform_glycosylation_sites_literature_mining.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))


    seen = {}
    file_list = glob.glob("downloads/lit_min/????_??_??/glycosylation_sites_*.csv")
    

    for in_file in file_list:
        file_name = in_file.split("/")[-1].replace(".csv", "")
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, ",")
        f_list = data_frame["fields"]
        for row in data_frame["data"]:
            uniprotkb_ac = row[f_list.index("uniprotkb_ac")].split("-")[0]
            aa_pos = row[f_list.index("glycosylation_site")]
            amino_acid = row[f_list.index("amino_acid")]
            evidence = row[f_list.index("evidence")]
            llm_flag = ""
            if "llm_verification_flag" in f_list:
                llm_flag = row[f_list.index("llm_verification_flag")]
            if llm_flag == "llm_no":
                continue
            combo = "%s|%s|%s|%s" % (uniprotkb_ac, aa_pos, amino_acid, evidence)
            if combo not in seen:
                seen[combo] = {}
            mining_tool = file_name.replace("glycosylation_sites_","")
            seen[combo][mining_tool] = True

    seen_row = {}
    row_list = []
    file_list = glob.glob("downloads/lit_min/current/glycosylation_sites_*.csv")
    for combo in seen:
        uniprotkb_ac, aa_pos, amino_acid, evidence = combo.split("|")
        tool_list = list(seen[combo].keys())
        #if "GlycoSiteMiner" in tool_list:
        #    tool_list.remove("GlycoSiteMiner")
        if len(tool_list) == 0:
            continue
        mining_tools =  ";".join(tool_list) 
    
        aa_three = aa_format_dict["three"][amino_acid]
        glycosylation_type = aa_format_dict["glytype"][aa_three].lower()
        saccharide = ""
        if uniprotkb_ac not in ac2canon:
            continue
        canon = ac2canon[uniprotkb_ac]
        aa_one = aa_format_dict["one"][aa_three]
        ds = "%s_proteoform_glycosylation_sites_literature_mining" % (species)
        src_xref_key = "protein_xref_glygen_ds"
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""

        pair_list = [["protein_xref_pubmed", evidence]]
        for pair in pair_list:
            flag_list = []
            if glycosylation_type == "":
                flag_list.append("bad_glytype")
                g_type = ""
            else:
                g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()

            #combo_id = "%s|%s" % (canon, aa_pos)
            #if combo_id not in known_site_dict:
            #    flag_list.append("site_not_reported")
            newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                pair[0],pair[1], src_xref_key, src_xref_id ,mining_tools]
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
                seen_row[newrowstr] = True
    FL.close()

    return out_df


               

def extract_glycosylation_sites_literature_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()

    
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
            if f not in required_output_fields and f not in extra_fields:
                if f not in ["cell_line_name", "cell_line_species"]:
                    extra_fields.append(f)


    out_df = {"fields":[], "data":[]}
    headers = required_output_fields + extra_fields
    newrow = headers
    out_df["fields"] = newrow


    FL1 = open("logs/%s_proteoform_glycosylation_sites_literature.1.log" % (species), "w")
    FL1.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    FL2 = open("logs/%s_proteoform_glycosylation_sites_literature.2.log" % (species), "w") 

    ds = "%s_proteoform_glycosylation_sites_literature" % (species)
    src_xref_key = "protein_xref_data_submission" if species == "hcv1a" else "protein_xref_glygen_ds"
    src_xref_id =  ds2bco[ds] if ds in ds2bco else ""

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
            xref_id = xref_id.replace("GLYDS", "GLY_")

            if xref_id.find("GLY_") != -1:
                xref_key = "protein_xref_unicarbkb_ds"


            aa_three = row[f_list.index("amino_acid")]
            glycosylation_type = row[f_list.index("glycosylation_type")]
            saccharide, cell_line, cellosaurus_id, abundance_normalized = "", "", "", ""
           
            if "saccharide" in f_list:
                saccharide = row[f_list.index("saccharide")]
            
            extra_values = []
            for f in extra_fields:
                if f not in ["cell_line_name", "cell_line_species"]:
                    v = row[f_list.index(f)] if f in f_list else ""
                    if f == "source_cell_line_cellosaurus_name":
                        if "source_cell_line_cellosaurus_id" in f_list:
                            clid = row[f_list.index("source_cell_line_cellosaurus_id")]
                            v = clid2name[clid] if clid in clid2name else ""
                    #if f == "source_cell_line_cellosaurus_id":
                    #    v = v.replace("_", ":")
                    extra_values.append(v)

            #RRRRRRR

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
                        if species in ["hcv1a"]:
                            newrow = [ac.split("-")[0], pos, aa_three,saccharide,g_type,
                                "protein_xref_unicarbkb_ds", src_xref_id, src_xref_key,src_xref_id] + extra_values
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
        
        #If saccharide is not in glycan_dict, force it to be ""
        saccharide = newrow[3].strip()
        if saccharide not in logged_glycan and saccharide not in glycan_dict:
            logged_glycan[saccharide] = True
            FL2.write("\"%s\"\n" % (saccharide))
            newrow[3] = ""
        out_df["data"].append(newrow)

    FL1.close()
    FL2.close()



     #take care of rows with unmapped tissueid and clid 
    log_file = "logs/%s_proteoform_glycosylation_sites_literature.3.log" % (species)
    FW = open(log_file, "w")
    FW.write("\"%s\"\n" % ("\",\"".join(out_df["fields"] + ["flags"])))
    row_list = []
    for newrow in out_df["data"]:
        # log unmappable tissueid 
        flag_list = []
        if "source_tissue_id" in out_df["fields"]:
            tid = newrow[out_df["fields"].index("source_tissue_id")]
            tid = "" if tid.strip() == "0" else tid
            tid = tid.replace("UBERON_", "UBERON:")
            tname = tissueid2name[tid] if tid in tissueid2name else ""
            newrow[out_df["fields"].index("source_tissue_id")] = tid
            newrow[out_df["fields"].index("source_tissue_name")] = tname
            if tid.strip() != "" and tname == "":
                flag_list.append("unmapped_source_tissue_id(%s)" %(tid))

        # log unmappable clid
        if "source_cell_line_cellosaurus_id" in out_df["fields"]:
            clid = newrow[out_df["fields"].index("source_cell_line_cellosaurus_id")]
            clid = "" if clid.strip() == "0" else clid
            cname = clid2name[clid] if clid in clid2name else ""
            newrow[out_df["fields"].index("source_cell_line_cellosaurus_name")] = cname
            if clid.strip() != "" and cname == "":
                flag_list.append("unmapped_source_cell_line_cellosaurus_id(%s)"%(clid))

        if flag_list != []:
            FW.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
        else:
            row_list.append(newrow)
    FW.close()
    out_df["data"] = row_list

    return out_df


def extract_glycosylation_sites_o_gluc_ds(species, ds_name):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()


    extra_fields_one = ["glycosylation_subtype", "carb_name", "notes"]
    newrow = []
    for v in required_output_fields + extra_fields_one:
        v = v.strip().replace(" ", "_").replace("(", "_").replace(")", "_").lower()
        v = v.replace("__", "_").replace(",","").replace("-","_").replace("/","_")
        if v[-1] == "_":
            v = v[:-1]
        newrow.append(v)

    out_df = {"fields":[], "data":[]}
    out_df["fields"] = newrow
    data_frame = {}
    in_file = "compiled/%s_proteoform_glycosylation_sites_o-gluc.csv"  % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    #['uniprot_kb_ac', 'gene_name', 'glycosylation_site', 'amino_acid', 'peptide_start', 'peptide_end', 'peptide_seq', 'saccharide', 'glycosylation_type', 'glycosylation_subtype', 'glycosylation stoichiometry', 'evidence', 'data_source', 'notes']
    
        
    FL = open("logs/%s_proteoform_glycosylation_sites_o_gluc.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    for row in data_frame["data"]:
        flag_list = []
        ac = row[f_list.index("uniprot_kb_ac")] 
        pos = row[f_list.index("glycosylation_site")] 
        evidence = row[f_list.index("evidence")]
        data_src = row[f_list.index("data_source")]
        notes = row[f_list.index("notes")]
        notes += " " + row[f_list.index("glycosylation_stoichiometry")]
        if ds_name == "glycosylation_sites_o_gluc" and data_src.lower() != "pubmed":
            continue
        if ds_name == "glycosylation_sites_o_gluc_predicted" and data_src.lower() != "prosite":
            continue
        if ds_name == "glycosylation_sites_o_gluc_predicted":
            notes = "Predicted O-glucosylation site, and attached O-Gluc glycan (G71142DF), based upon the presence of the consensus sequence C3-x-N-T-x-G-S-(FY)-x-C4 for the O-glucosyltransferases POGLUT2 and POGLUT3 identified in their paper (PMID:34411563)."

        #protein_xref_prosite-prorule, protein_xref_glygen_ds
        canon = ac2canon[ac] if ac in ac2canon else ""
        src_xref_key = "protein_xref_data_submission" if species == "human" else "protein_xref_glygen_ds" 
        ds = "%s_proteoform_glycosylation_sites_o_gluc" % (species)
        src_xref_id = ds2bco[ds] if ds in ds2bco else ""
        tmp_dict = {
            "uniprotkb_canonical_ac":canon,
            "src_xref_key":src_xref_key, "src_xref_id":src_xref_id,
            "glycosylation_site_uniprotkb":pos,
            "notes":notes
        }
        if data_src.lower() == "prosite":
            tmp_dict["xref_key"] = "protein_xref_prosite_prorule"
            tmp_dict["xref_id"] = "PRU00076"
        elif data_src.lower() == "pubmed":
            if evidence.strip() != "":
                tmp_dict["xref_key"] = "protein_xref_pubmed"
                tmp_dict["xref_id"] = evidence
        
        flag_list = []
        if ac not in ac2canon:
            flag_list.append("not mapping to canonincal")
        val_dict = {}
        newrow = []
        for f in required_output_fields + extra_fields_one:
            v = row[f_list.index(f)] if f in f_list else "xxxxx"
            v = tmp_dict[f] if f in tmp_dict else v
            val_dict[f] = v
            newrow.append(v)
        
        aa_three = val_dict["amino_acid"] 
        g_type = val_dict["glycosylation_type"]
        flag_list += get_gly_filter_flags(species,canon,pos,aa_three, g_type)
        if flag_list != []:
            FL.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
        else:
            out_df["data"].append(newrow)
            if ds_name == "glycosylation_sites_o_gluc_predicted":
                extrarow = json.loads(json.dumps(newrow))
                extrarow[out_df["fields"].index("xref_key")] = "protein_xref_pubmed"
                extrarow[out_df["fields"].index("xref_id")] = "34411563"
                out_df["data"].append(extrarow)

    FL.close()


    return out_df


def extract_glycosylation_sites_oglcnac_mcw_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()



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

    
    data_frame = {}
    in_file = "downloads/mcw_oglcnac/current/%s_o_glcnacome_mcw.csv" % (species)
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
    for v in required_output_fields + extra_fields_one + extra_fields_two:
        v = v.strip().replace(" ", "_").replace("(", "_").replace(")", "_").lower()
        v = v.replace("__", "_").replace(",","").replace("-","_").replace("/","_")
        if v[-1] == "_":
            v = v[:-1]
        newrow.append(v)

    out_df = {"fields":[], "data":[]}
    out_df["fields"] = newrow
    FL = open("logs/%s_proteoform_glycosylation_sites_oglcnac_mcw.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))


    #ac_field = "ENTRY ID"
    ac_field = "UniprotKB ID"
    #sites_field = "O-GlcNAc Sites"
    sites_field = "oglcnac sites"
    score_field = "oglcnacscore"
    pmid_field = "PMIDS"

    row_list = {}
    seen = {}
    seen_row = {}
    for row in data_frame["data"]:
        flag_list = []
        ac = row[f_list.index(ac_field)]
        if row[f_list.index(sites_field)].strip() == "":
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
        pmid = row[f_list.index(pmid_field)].strip()
        if pmid != "":
            pmid_list += pmid.strip().split(";")
        pmid_list = list(set(pmid_list))

        oglcnac_score = row[f_list.index(score_field)].strip()
        site_str = row[f_list.index(sites_field)].strip()
        site_str = site_str.replace("(", "").replace(")", "")
        site_str = site_str.replace(" or ", " / ")
        site_list = site_str.split(";")
        aa_pos = ""
        aa_three = ""
        glycosylation_type = "o-linked"
        saccharide = "G49108TO"
        if saccharide not in glycan_dict:
            flag_list.append("saccharide not in glycan list")


        for s in site_list:
            aa_one, aa_three, aa_pos = "","", ""
            if s.strip() != "":
                s = s.strip()
                aa_one, aa_pos = s[0], s[1:].split(" ")[0]
                aa_three = aa_format_dict["three"][aa_one]

            src_xref_key = "protein_xref_oglcnac_db"
            src_xref_id = canon.split("-")[0]
            for pmid in pmid_list:
                pair_list = [["protein_xref_pubmed", pmid]]
                for pair in pair_list:
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [canon, str(aa_pos), aa_three, saccharide, g_type,
                            pair[0],pair[1], src_xref_key, src_xref_id] 
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


    return out_df



def extract_glycosylation_sites_oglcnac_atlas_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()

    in_file = "generated/misc/sample_mapping.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    #term_id,term_name,sample_type,ontology_name
    sampletype2ontology = {}
    for row in data_frame["data"]:
        term_id,term_name = row[f_list.index("term_id")],  row[f_list.index("term_name")]
        sample_type_list = [row[f_list.index("term_name")].strip()]
        #for v in row[f_list.index("in_dataset")].split(","):
        in_dataset = row[f_list.index("in_dataset")].strip()
        sample_type_list.append(in_dataset)
        for sample_type in sample_type_list:
            sample_type = sample_type.lower()
            if sample_type not in sampletype2ontology:
                sampletype2ontology[sample_type] = {}
            sampletype2ontology[sample_type][term_id] = term_name
    

    extra_fields_two = ["carb_name", "glycosylation_subtype"]
    extra_fields_three = [
        "source_tissue_id", "source_tissue_name",
        "source_cell_line_cellosaurus_id", "source_cell_line_cellosaurus_name"
    ]        
    #extra_fields_four = ["ontology_term_id","ontology_term_name"]
    
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
            #extra_dict_two[canon]["eco_id"] = "ECO_0000269"
            extra_dict_two[canon]["carb_name"] = "GlcNAc"
            extra_dict_two[canon]["glycosylation_subtype"] = "O-GlcNAcylation"

            for f in f_list:
                if f in extra_fields_two and f not in extra_dict_two[canon]:
                    extra_dict_two[canon][f] = row[f_list.index(f)]

    
    data_frame = {}
    in_file_one = "downloads/atlas_oglcnac/current/unambiguous_sites.csv"
    in_file_two = "downloads/atlas_oglcnac/current/ambiguous_sites.csv"
    #,id,species,sample_type,accession,accession_source,entry_name,protein_name,gene_name,peptide_seq,site_residue,position_in_peptide,position_in_protein,method,analytical_throughput,pmid,comments
    
    data_frame, data_frame_one, data_frame_two = {}, {}, {}
    libgly.load_sheet(data_frame_one, in_file_one, ",")
    libgly.load_sheet(data_frame_two, in_file_two, ",")



    data_frame["fields"] = data_frame_one["fields"]
    data_frame["data"] = data_frame_one["data"] + data_frame_two["data"]
    #data_frame["data"] = data_frame_one["data"]

    

    ignored_fields = ["entry_name","protein_name","id","species","status","gene_name","pmid"]
    ignored_fields += ["", "accession", "accession_source", "peptide_seq"]
    ignored_fields += ["site_residue","position_in_peptide", "position_in_protein" ]
    
    
    f_list = data_frame["fields"]
    extra_fields_one = []
    for f in f_list:
        if f.lower() in ignored_fields:
            continue
        f = "notes" if f == "comments" else f
        extra_fields_one.append(f)

    newrow = []
    for v in required_output_fields + extra_fields_one + extra_fields_two:
        newrow.append(v)
    #for v in extra_fields_three + extra_fields_four:
    for v in extra_fields_three:
        newrow.append(v)
    newrow.append("peptide")

    out_df = {"fields":[], "data":[]}
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_glycosylation_sites_oglcnac_atlas.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))


    #ac_field = "ENTRY ID"
    ac_field = "accession"
    #sites_field = "O-GlcNAc Sites"
    pos_field = "position_in_protein"
    aa_field = "site_residue"
    pmid_field = "pmid"

    row_list = {}
    seen = {}
    seen_row = {}
    for row in data_frame["data"]:
        flag_list = []
        ac = row[f_list.index(ac_field)]
        #sample_type = row[f_list.index("sample_type")].replace("cells", "").strip()
        sample_type = row[f_list.index("sample_type")].strip()
        canon = ac2canon[ac] if ac in ac2canon else ""
        if row[f_list.index(pos_field)].strip() == "":
            flag_list.append("no o-site position")
        if row[f_list.index(aa_field)].strip() == "":
            flag_list.append("no o-site amino acid")
        #if ac not in ac2canon:
        #    flag_list.append("not mapping to canonincal")
        
        extra_values_one = []
        for f in extra_fields_one:
            f = "comments" if f == "notes" else f
            extra_values_one.append(row[f_list.index(f)].replace("\n", ""))
        
        extra_values_two = []
        for f in extra_fields_two:
            v = ""
            if canon in extra_dict_two:
                v = extra_dict_two[canon][f] if f in extra_dict_two[canon] else ""
            v = v.replace("\n", "")
            extra_values_two.append(v)

        pmid_list = []
        pmid = row[f_list.index(pmid_field)].strip()
        if pmid != "":
            pmid_list += pmid.strip().split(";")
        pmid_list = list(set(pmid_list))
        aa_pos = row[f_list.index(pos_field)].strip()
        aa_one = row[f_list.index(aa_field)].strip()
        aa_three = aa_format_dict["three"][aa_one] if aa_one in aa_format_dict["three"] else ""

        glycosylation_type = "o-linked"
        glycosylation_subtype = "O-GlcNAcylation"
        saccharide = "G49108TO"
        if saccharide not in glycan_dict:
            flag_list.append("saccharide not in glycan list")
        if aa_pos.isdigit() == False:
            flag_list.append("bad_pos=%s" %(aa_pos))

        src_xref_key = "protein_xref_oglcnac_atlas"
        src_xref_id = canon.split("-")[0]
        for pmid in pmid_list:
            pair_list = [["protein_xref_pubmed", pmid]]
            if pmid.find("doi") != -1:
                pmid = pmid.replace("doi.org/","")
                pair_list = [["protein_xref_doi", pmid]]
            for pair in pair_list:
                g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                newrow = [canon, str(aa_pos), aa_three, saccharide, g_type,
                            pair[0],pair[1], src_xref_key, src_xref_id] 
                newrow += extra_values_one + extra_values_two
                
                #relax site requirement
                if "no o-site" in flag_list:
                    flag_list.remove("no o-site")

                if flag_list  != []:
                    newrow[0] = ac
                    if ac in ac2canon:
                        #output species specific rows in logfile
                        FL.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                else:
                    flag_list = get_gly_filter_flags(species,canon,aa_pos,aa_three, glycosylation_type)
                    #relax no-site rows
                    if aa_pos == "" and aa_three == "" and "aa_not_in_dictionary" in flag_list:
                        flag_list.remove("aa_not_in_dictionary")
                    if flag_list  != []:
                        #output species specific rows in logfile
                        if ac in ac2canon:
                            FL.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                    else:
                        peptide = ""
                        if canon != "" and aa_pos != "":
                            aa_idx = int(aa_pos)-1
                            pep_len = 5
                            pep_one = seq_hash[canon][:aa_idx][-pep_len:]
                            pep = seq_hash[canon][aa_idx]
                            pep_two = seq_hash[canon][int(aa_pos):][:pep_len]
                            peptide = pep_one +"-"+ pep + "-" + pep_two
                                          
                        if canon not in row_list:
                            row_list[canon] = {}
                        if aa_pos not in row_list[canon]:
                            row_list[canon][aa_pos] = []
                        onto_list = []
                        sample_type_lc = sample_type.lower()
                        if sample_type_lc in sampletype2ontology:
                            for term_id in sampletype2ontology[sample_type_lc]:
                                term_name = sampletype2ontology[sample_type_lc][term_id]
                                onto_list.append([term_id, term_name])
                        onto_list = [["", ""]] if onto_list == [] else onto_list
                        for r in onto_list:
                            #r[0] = r[0].replace("_", ":")
                            rr = r + ["", ""] 
                            if r[0].find("CVCL_") != -1:
                                rr = ["",""] + r
                            row_list[canon][aa_pos].append(newrow + rr + [peptide])
    FL.close()


    seen = {}
    for canon in sorted(row_list.keys()):
        for aa_pos in sorted(row_list[canon].keys()):
            for newrow in row_list[canon][aa_pos]:
                row_str = json.dumps(newrow)
                if row_str not in seen:
                    out_df["data"].append(newrow)
                    seen[row_str] = True


    return out_df


def extract_glycosylation_sites_tyr_o_linked_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()


    data_frame = {}
    in_file = "compiled/%s_proteoform_glycosylation_sites_tyr_o_linked.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]

    out_df = {"fields":[], "data":[]}
    newrow = required_output_fields
    out_df["fields"] = newrow
    FL = open("logs/%s_proteoform_glycosylation_sites_tyr_o_linked.log" % (species), "w")
    FL.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    seen = {}
    seen_row = {}
    for row in data_frame["data"]:
        ac = row[f_list.index("uniprotkb_ac")]
        aa_pos = row[f_list.index("glycosylation_site_uniprotkb")]
        aa_three = row[f_list.index("amino_acid_uniprotkb")]
        glycosylation_type = row[f_list.index("glycosylation_type")].lower()
        evidence = row[f_list.index("evidence")].replace("GLYDS", "GLY_")

        pmid = row[f_list.index("pmid")]
        if ac not in ac2canon:
            continue
        canon = ac2canon[ac] if ac in ac2canon else ac
        saccharide = ""
        pair_list = [["protein_xref_pubmed", pmid]]
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
                seen_row[newrowstr] = True

    FL.close()

    return out_df





def extract_glycosylation_sites_harvard_ds(species):

    
    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()

    data_frame = {}
    in_file = path_obj["downloads"] + "harvard/%s_glycosylation_current.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
   
    out_df = {"fields":[], "data":[]}
    newrow = required_output_fields + ["composition"]
    out_df["fields"] = newrow
    FL1 = open("logs/%s_proteoform_glycosylation_sites_harvard.1.log" % (species), "w")
    FL1.write("%s\n" % ("\",\"".join(newrow + ["filter_flags"])))

    FL2 = open("logs/%s_proteoform_glycosylation_sites_harvard.2.log" % (species), "w")


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
                evidence = evidence.replace("GLYDS", "GLY_")
                pair_list = [["protein_xref_pubmed", pmid]]
                for pair in pair_list:
                    #If saccharide is not in glycan_dict, force it to be ""
                    if saccharide not in logged_glycan and saccharide.strip() != "" and saccharide not in glycan_dict:
                        logged_glycan[saccharide] = True
                        FL2.write("\"%s\"\n" % (saccharide))

                        saccharide = ""
                    src_xref_key = "protein_xref_data_submission" if species == "human" else "protein_xref_glygen_ds"
                    g_type = glycosylation_type[0].upper() + glycosylation_type[1:].lower()
                    newrow = [canon,str(aa_pos),aa_three,saccharide,g_type,
                            pair[0],pair[1],src_xref_key,evidence, compo]
                    if flag_list != []:
                        FL1.write("%s\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
                        continue
                    newrowstr = json.dumps(newrow)
                    if newrowstr not in seen:
                        out_df["data"].append(newrow)
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




def load_glytoucan_type_dict(gtc_type_dict):


    data_frame = {}
    in_file = "unreviewed/glycan_masterlist.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    tmp_dict = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")].strip()
        glytoucan_type = row[f_list.index("glytoucan_type")].strip().lower()
        if glytoucan_ac not in tmp_dict:
            tmp_dict[glytoucan_ac] = {}
        tmp_dict[glytoucan_ac][glytoucan_type] = True


    for glytoucan_ac in tmp_dict:
        gtc_types = ";".join(sorted(list(tmp_dict[glytoucan_ac].keys())))
        gtc_type_dict[glytoucan_ac] = gtc_types


    return

        
        

def load_glycosylation_type_two():
 
     
    data_frame = {}
    in_file = "unreviewed/glycan_classification.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    glytoucanac2glycosylationtype = {}
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")].strip()
        gly_type = row[f_list.index("glycan_type")].strip().lower()
        if glytoucan_ac not in glytoucanac2glycosylationtype:
            glytoucanac2glycosylationtype[glytoucan_ac] = []
        if gly_type not in glytoucanac2glycosylationtype[glytoucan_ac]:
            if gly_type == "n-linked":
                glytoucanac2glycosylationtype[glytoucan_ac].append(gly_type)
            if gly_type == "o-linked":
                glytoucanac2glycosylationtype[glytoucan_ac].append(gly_type)

    data_frame = {}
    in_file = "generated/misc/glycan_add_class.csv"
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("glytoucan_ac")].strip()
        gly_type = row[f_list.index("glycan_type")].strip().lower()
        if glytoucan_ac not in glytoucanac2glycosylationtype:
            glytoucanac2glycosylationtype[glytoucan_ac] = []
        if gly_type not in glytoucanac2glycosylationtype[glytoucan_ac]:
            if gly_type in ["n-linked", "o-linked", "c-linked"]:
                glytoucanac2glycosylationtype[glytoucan_ac].append(gly_type)

    
 
    return glytoucanac2glycosylationtype

def load_glycan_dict():
    glycan_dict = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/glycan_properties.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        ac = row[f_list.index("glytoucan_acc")]
        glycan_dict[ac] = True

    return glycan_dict






def map_disease_info (obj, do_id):

    obj["do_id"] = do_id
    if obj["do_id"].find(":") != -1:
        obj["do_id"] = obj["do_id"].split(":")[1]
    if obj["do_id"] in doid2name:
        obj["disease_name"] = doid2name[obj["do_id"]]

    return


def map_tissue_info (obj, uberon_id):

    obj["uberon_id"] = uberon_id
    if obj["uberon_id"].find(":") != -1:
        obj["uberon_id"] = obj["uberon_id"].split(":")[1]
    if obj["uberon_id"] in uberonid2name:
        obj["tissue_name"] = uberonid2name[obj["uberon_id"]]

    return


def load_disease_info ():

    dict_one, dict_two = {}, {}
    in_file = "generated/misc/doid2uberonid_mapping.csv"
    
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        dict_one[row[f_list.index("do_id")]] = row[f_list.index("do_name")]
        dict_two[row[f_list.index("uberon_id")]] = row[f_list.index("uberon_name")]

    return dict_one, dict_two


def load_cellosaurus_info():

    in_file = "downloads/cellosaurus/current/cellosaurus.txt"
    log_non_csv_file_use(in_file)

    cellosaurus_info = {}
    with open(in_file, "r") as FR:
        for line in FR:
            if line[0:2] == "ID":
                cl_name = line[3:].strip()
            if line[0:2] == "AC":
                cl_id = line[3:].strip()
                cellosaurus_info[cl_id] = cl_name
  
    append_from_sample_mapping(cellosaurus_info)
    return cellosaurus_info


def append_from_sample_mapping(in_dict):
    
    in_file = "generated/misc/sample_mapping.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        in_src = row[f_list.index("in_dataset")]
        term_name = row[f_list.index("term_name")]
        in_dict[in_src] = term_name
    return 



def  load_uckb2glytoucan():

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

    return uckb2glytoucan


def load_carbbank2glytoucan():

    tmp_dict = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/carbbank.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        glytoucan_ac = row[f_list.index("GlyTouCanAccession")]
        tmp_id = row[f_list.index("CarbbankAccession")]
        if tmp_id not in tmp_dict:
            tmp_dict[tmp_id] = []
        tmp_dict[tmp_id].append(glytoucan_ac)

    return tmp_dict



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
    append_from_sample_mapping(tmp_dict)

    return tmp_dict


def load_clid2name_map():

    in_file = "downloads/cellosaurus/current/cellosaurus.txt"
    log_non_csv_file_use(in_file)
    tmp_dict = {}
    with open(in_file, "r") as FR:
        for line in FR:
            if line[0:2] == "ID":
                cl_name = line[3:].strip()
            if line[0:2] == "AC":
                cl_id = line[3:].strip()
            if line[0:2] == "//":
                tmp_dict[cl_id] = cl_name
                cl_id = cl_id.replace("_", ":")
                tmp_dict[cl_id] = cl_name


    in_file = "generated/misc/sample_mapping.csv"
    data_frame = {}
    libgly.load_sheet(data_frame, in_file, ",")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        term_id,term_name = row[f_list.index("term_id")],  row[f_list.index("term_name")]
        sample_type = row[f_list.index("sample_type")]
        if sample_type == "cell_line":
            cl_id = term_id
            tmp_dict[cl_id] = term_name
            cl_id = cl_id.replace("_", ":")
            tmp_dict[cl_id] = term_name

    append_from_sample_mapping(tmp_dict)
    return tmp_dict





def extract_glycosylation_sites_unicarbkb_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()

    exclude_dict = {}
    line_list = open("generated/misc/exclude_list.csv", "r").read().split("\n")
    for line in line_list:
        if line.strip() != "":
            ac, pmid = line.split(",")[0].replace("\"", ""), line.split(",")[-1].replace("\"", "")
            combo_id = "%s|%s" % (ac, pmid)
            exclude_dict[combo_id] = True


    required_input_fields = ["protein","position", "aminoacid", "toucan","aminoacidtype"]
    append_fields = ["saccharide"]
    extra_fields = [
        "uckb_id","composition","curation_notes","additional_notes",
        "pdb", "swiss_model", "abundance",
        "source_tissue_id", "source_tissue_name",  
        "source_cell_line_cellosaurus_id", "source_cell_line_cellosaurus_name"
    ]

    f_map = { 
        "acc":"protein",
        "glycosylation_type":"aminoacidtype",
        "site":"position",
        "protein":"uniprotkb_ac",
        "id":"uckb_id",
        "toucan":"saccharide",
        "glytoucan":"saccharide",
        "glycosite":"position", 
        "amino_acid":"aminoacid",
        "typeaminoacid":"glycotype", 
        "type_amino_acid":"glycotype",
        "type":"glycotype",
        "aminoacidtype":"glycotype",
        "additionalnotes":"additional_notes",
        "notes":"curation_notes",
        "pos_start":"start_pos",
        "pos_end":"end_pos",
        "peptide_sequence":"peptide",
        "peptide_start":"start_pos",
        "peptide_end":"end_pos"
    }


    out_df = {"fields":[], "data":[]}
    out_df["fields"] = required_output_fields
    out_df["fields"] += ["start_pos","end_pos", "start_aa","end_aa","site_seq", "site_type", "src_file_name"]
    # append mapped field names
    out_df["fields"] += extra_fields

    cellosaurus_info = load_cellosaurus_info()
    #print json.dumps(cellosaurus_info, indent=4)
    #exit()

    file_list = glob.glob(path_obj["downloads"] + "unicarbkb/current/*.csv")
    file_list += glob.glob(path_obj["downloads"] + "unicarbkb/current/*.tsv")



    input_obj_list = []
    log_dict = {}
    for in_file in file_list:
        file_name = ".".join(in_file.split("/")[-1].split(".")[:-1])
        file_ext = in_file.split(".")[-1].lower()
        sep = "\t" if file_ext == "tsv" else ","
        data_frame = {}
        libgly.load_sheet(data_frame, in_file, sep)
        f_list = data_frame["fields"]
        for i in xrange(0,len(f_list)):
            f_list[i] = f_list[i].lower()
            if f_list[i]  in ["tissue", "tissue_id"]:
                f_list[i] = "source_tissue_id"
            if f_list[i]  in ["cellosaurus"]:
                f_list[i] = "source_cell_line_cellosaurus_id"
            for f in f_map:
                f_list[i] = f_map[f] if f_list[i] == f else f_list[i]
        
        #append some fields to pass required_input_fields criteria
        for fa in append_fields:
            if fa not in f_list:
                f_list.append(fa)

        
        for f in required_input_fields:
            if f not in f_list and f_map[f] not in f_list:
                print "Required field %s not in f_list for %s" % (f.lower(), file_name)
                print f_list
                exit()
        row_idx = 0
        in_species = 0
   


        for row in data_frame["data"]:
            #append values for fa in append_fields
            for fa in append_fields:
                row.append("")
            row_idx += 1
            ac = row[f_list.index("uniprotkb_ac")]
            sp = ""
            if ac.strip() == "":
                if species_obj[species]["long_name"].lower() == row[f_list.index("species")].lower():
                    sp = species
            elif ac in ac2canon_strict:
                sp = species
            if sp == species:
                in_species += 1
                obj = {"filename":file_name, "rowidx":row_idx, "origrow":row}
                for f in f_list:
                    obj[f] = row[f_list.index(f)].split("^")[0] if f == "position" else row[f_list.index(f)]
                    obj[f] = obj[f].strip()
                

                if obj["aminoacid"].strip() != "": 
                    obj["aminoacid"] = obj["aminoacid"].split("_")[-1]
                elif obj["glycotype"].lower() in ["asn", "thr", "ser"]:
                    obj["aminoacid"] = obj["glycotype"]

                for f in ["composition", "curation_notes", "additional_notes"]:
                    obj[f] = row[f_list.index(f)].strip() if f in f_list else ""

                input_obj_list.append(obj)

                obj["do_id"], obj["uberon_id"], obj["disease_name"], obj["tissue_name"], = "","","",""
                if "doid" in f_list:
                    map_disease_info(obj, row[f_list.index("doid")].strip())
                if "uberon_id" in f_list:
                    map_tissue_info(obj, row[f_list.index("uberon_id")].strip())
            else:
                if file_name not in log_dict:
                    log_row = ["row_num"] + ["flag_list", "qc_call"]
                    log_dict[file_name] = [ log_row]
                #log_row = [str(row_idx)] + row + ["not_in_species","qc_failed"]
                log_row = [str(row_idx)] + ["not_in_species","qc_failed"]
                log_dict[file_name].append(log_row)


    out_obj_list = []
    for in_obj in input_obj_list:
        if "tissue_id" in in_obj:
            w_list = in_obj["tissue_id"].strip().split(" ")
            if len(w_list) > 1:
                for w in w_list:
                    if w.find("UBERON:") != -1:
                        in_obj["tissue_id"] = w.replace(";", "")
                        break

        in_obj["evdnlist"] = []
        if "pmid" in in_obj:
            if in_obj["pmid"] not in ["", "0"]:
                o = {"xref_key":"protein_xref_pubmed", "xref_id":in_obj["pmid"]}
                in_obj["evdnlist"].append(o)
       
        if "doi" in in_obj:
            if in_obj["doi"] != "":
                in_obj["doi"] = in_obj["doi"].replace("https://doi.org/", "")
                in_obj["doi"] = in_obj["doi"].replace("doi.org/", "")
                o = {"xref_key":"protein_xref_doi", "xref_id":in_obj["doi"]}
                in_obj["evdnlist"].append(o)
       
        in_obj["gtclist"] = [] if in_obj["saccharide"] == "" else in_obj["saccharide"].split("|")
        if "uckb_id" in in_obj:
            if in_obj["saccharide"] == "" and in_obj["uckb_id"][0:5] == "comp_":
                in_obj["composition"] = in_obj["uckb_id"]
                if in_obj["composition"] in uckb2glytoucan:
                    in_obj["gtclist"] += uckb2glytoucan[in_obj["composition"]]

            if in_obj["uckb_id"][0:5] == "comp_":
                in_obj["uckb_id"] = ""



        in_obj["cl_name"], in_obj["clo_id"], in_obj["clo_name"] = "", "", ""
        if "cellosaurus" in in_obj:
            cl_id = in_obj["cellosaurus"]
            if cl_id in cellosaurus_info:
                in_obj["cl_name"] = cellosaurus_info[cl_id]


        in_obj["aalist"], in_obj["poslist"] = [], []
        if in_obj["aminoacid"].find("|") != -1:
            in_obj["aalist"] = in_obj["aminoacid"].split("|")
            in_obj["poslist"] = in_obj["position"].split("|")
        elif in_obj["position"].isdigit() == False and in_obj["position"].strip() != "":
            if in_obj["position"].find(";") != -1:
                site_list = in_obj["position"].split(";")
                for site in site_list:
                    if site[1:].isdigit() == True:
                        in_obj["aalist"].append(site[0])
                        in_obj["poslist"].append(site[1:])
            elif in_obj["position"][1:].isdigit() == True:
                if in_obj["position"][1:].isdigit() == True:
                    in_obj["aalist"] = [in_obj["position"][0]]
                    in_obj["poslist"] = [in_obj["position"][1:]]
            elif in_obj["position"].find("|") != -1:
                for site in in_obj["position"].split("|"):
                    if site.isdigit() == True:
                        in_obj["poslist"].append(site)
        elif in_obj["position"].isdigit() == True:
            in_obj["aalist"].append(in_obj["aminoacid"])
            in_obj["poslist"].append(in_obj["position"])


        #Process peptide range based 
        site_type = ""
        if in_obj["poslist"] != []:
            site_type = "known_site"
            in_obj["start_pos"] = int(in_obj["poslist"][0])
            in_obj["end_pos"] = int(in_obj["poslist"][0])
        elif "start_pos" in in_obj:
            in_obj["start_pos"] = in_obj["start_pos"].split(".")[0]
            in_obj["end_pos"] = in_obj["end_pos"].split(".")[0]
            if in_obj["start_pos"].isdigit() and in_obj["end_pos"].isdigit():
                in_obj["start_pos"] = int(in_obj["start_pos"])
                in_obj["end_pos"] = int(in_obj["end_pos"])
                site_type = "known_site" if in_obj["start_pos"] == in_obj["end_pos"] else "fuzzy_site"
        elif "peptide" in in_obj:
            site_type = "unknown_site"
            in_obj["start_pos"] = -1
            in_obj["end_pos"] = -1

        if site_type == "":
            continue
        
        in_obj["site_seq"] = ""
        if in_obj["start_pos"] != -1:
            in_obj["site_seq"] = "extract_it"
        elif "peptide" in in_obj:
            in_obj["site_seq"] = in_obj["peptide"]

        #At this point, we know what type of site it is 

        flag_dict = {}
        if in_obj["uniprotkb_ac"] not in ac2canon_strict:
            flag_dict["uniprotkb_ac_not_in_masterlist"] = True
        else:
            in_obj["canon"] = ac2canon_strict[in_obj["uniprotkb_ac"]] 
            if site_type in ["fuzzy_site"]:
                if "start_pos" in in_obj and "end_pos" in in_obj and "peptide" in in_obj:
                    qc_glyco_fuzzy_position(in_obj["start_pos"], in_obj["end_pos"], 
                            in_obj["peptide"], in_obj["canon"], flag_dict)
            else:
                qc_glyco_position(in_obj["position"], in_obj["aminoacid"], in_obj["canon"], flag_dict) #unicarbkb_ds
         
        #if site_type in ["known_site", "mixed_site"]:
        in_obj["glycotype"] = qc_glyco_type(in_obj["aminoacid"],in_obj["saccharide"],flag_dict)


        if in_obj["pmid"] in ["", "0"] and in_obj["doi"] in ["", "0"]:
            flag_dict["no_publication_ref"] = True
        if in_obj["saccharide"] != "" and in_obj["saccharide"] not in glycan_dict:
            flag_dict["glytoucan_not_in_masterlist"] = True                                     
        if in_obj["glycotype"] == "":
            flag_dict["unknown_glycotype"] = True  

        #this filter is based on error found for 34106099
        tmp_id = "%s|%s" % (in_obj["uniprotkb_ac"].split("-")[0], in_obj["pmid"])
        if tmp_id in exclude_dict:
            flag_dict["rat2mouse_confusion_from_34106099"] = True

        flag_list = flag_dict.keys()
        qc_call = "qc_passed" if flag_list == [] else "qc_failed"
        file_name = in_obj["filename"]
        if file_name not in log_dict:
            log_row = ["row_num"] + ["flag_list", "qc_call"]
            log_dict[file_name] = [log_row]
        log_row = [str(in_obj["rowidx"])] + in_obj["origrow"] + [";".join(flag_list), qc_call]
        #log_row = [str(in_obj["rowidx"])] +  [";".join(flag_list), qc_call]
        log_dict[file_name].append(log_row)


        if flag_list == []:
            g_type = in_obj["glycotype"][0].upper() + in_obj["glycotype"][1:].lower()
            pos_list = "|".join(in_obj["poslist"])
            aa_list = "|".join(in_obj["aalist"])
            gtc_list = "|".join(in_obj["gtclist"])
            start_pos, end_pos = str(in_obj["start_pos"]), str(in_obj["end_pos"])
            
            site_seq = in_obj["site_seq"]
            if site_seq == "extract_it":
                s, e = int(start_pos), int(end_pos)
                site_seq = seq_hash[in_obj["canon"]][s-1:e]
            start_aa, end_aa = "", ""
            if site_seq.strip() != "":
                s_aa, e_aa = site_seq[0], site_seq[-1]
                start_aa = aa_format_dict["three"][s_aa] if s_aa in aa_format_dict["three"] else ""
                end_aa = aa_format_dict["three"][e_aa] if e_aa in aa_format_dict["three"] else ""


            ds = "%s_proteoform_glycosylation_sites_unicarbkb" % (species)
            src_xref_id = ds2bco[ds] if ds in ds2bco else ""
            src_xref_key = "protein_xref_unicarbkb_ds"
            #src_xref_key, src_xref_id = "protein_xref_unicarbkb_ds", in_obj["uniprotkb_ac"]
            for evdn_obj in in_obj["evdnlist"]:
                xref_key, xref_id = evdn_obj["xref_key"], evdn_obj["xref_id"]
                newrow = [in_obj["canon"],pos_list, aa_list, gtc_list,g_type,
                    xref_key,xref_id,src_xref_key, src_xref_id, start_pos, end_pos, 
                    start_aa, end_aa, site_seq,site_type, file_name]
                
                if "source_tissue_id" in in_obj:
                    for t in in_obj["source_tissue_id"].split(";"):
                        t = t.strip()
                        # this case is source_cell_line_cellosaurus_id not source_tissue_id
                        # thus correcting it here
                        if t in tissueid2name:
                            t = t.replace("_", ":")
                            if t.find(":") != -1:
                                if t.split(":")[1].isdigit():
                                    in_obj["source_tissue_id"] = t
                                    break
                    tmp_idx = in_obj["source_tissue_id"].find("UBERON:")
                    if tmp_idx != -1:
                        uberon_id = in_obj["source_tissue_id"][tmp_idx:].split(" ")[0].strip()
                        in_obj["source_tissue_id"] = uberon_id
                    if in_obj["source_tissue_id"] in tissueid2name:
                        in_obj["source_tissue_name"] = tissueid2name[in_obj["source_tissue_id"]]
                if "source_cell_line_cellosaurus_id" in in_obj:
                    #in_obj["source_cell_line_cellosaurus_id"] = in_obj["source_cell_line_cellosaurus_id"].replace("_", ":")
                    if in_obj["source_cell_line_cellosaurus_id"] in clid2name:
                        in_obj["source_cell_line_cellosaurus_name"] = clid2name[in_obj["source_cell_line_cellosaurus_id"]]

                for f in extra_fields:
                    v = in_obj[f] if f in in_obj else ""
                    newrow.append(v)
                out_df["data"].append(newrow)


    for file_name in log_dict:
        log_file = path_obj["logs"] + "unicarbkb/%s_%s.log" % (species,file_name)
        with open(log_file, "w") as FW:
            for row in log_dict[file_name]:
                FW.write("\"%s\"\n" % ("\",\"".join(row)))


    #take care of rows with unmapped tissueid and clid 
    log_file = "logs/%s_proteoform_glycosylation_sites_unicarbkb.3.log" % (species)
    FW = open(log_file, "w")
    FW.write("\"%s\"\n" % ("\",\"".join(out_df["fields"] + ["flags"])))
    row_list = []
    for newrow in out_df["data"]:
        flag_list = []
        # log unmappable tissueid 
        tid = newrow[out_df["fields"].index("source_tissue_id")]
        tid = "" if tid.strip() == "0" else tid
        tid = tid.replace("UBERON_", "UBERON:")
        if tid[:3] == "CL_":
            tid = tid.replace("CL_", "CL:")
        tname = tissueid2name[tid] if tid in tissueid2name else ""
        newrow[out_df["fields"].index("source_tissue_id")] = tid
        newrow[out_df["fields"].index("source_tissue_name")] = tname
        if tid.strip() != "" and tname == "":
            flag_list.append("unmapped_source_tissue_id(%s)" %(tid))
        
        # log unmappable clid
        clid = newrow[out_df["fields"].index("source_cell_line_cellosaurus_id")]
        clid = "" if clid.strip() == "0" else clid
        cname = clid2name[clid] if clid in clid2name else ""
        newrow[out_df["fields"].index("source_cell_line_cellosaurus_name")] = cname
        if clid.strip() != "" and cname == "":
            flag_list.append("unmapped_source_cell_line_cellosaurus_id(%s)"%(clid))
       
        if flag_list != []:
            FW.write("\"%s\"\n" % ("\",\"".join(newrow + [";".join(flag_list)])))
        else:
            row_list.append(newrow)
    FW.close()
    out_df["data"] = row_list

    return out_df




def qc_glyco_type(amino_acid,saccharide,flag_dict):

    aa_three = amino_acid[0].upper() +  amino_acid[1:].lower() if amino_acid != "" else ""
    glycosylation_type = "xxx: should change by under one of the conditions"
    if saccharide == "":
        if amino_acid.find("|") != -1:
            aa_list = amino_acid.split("|")
            f = aa_list[0]
            if aa_list[0] in aa_format_dict["glytype"]:
                glycosylation_type = aa_format_dict["glytype"][aa_list[0]]
            else:
                glycosylation_type = ""
                flag_dict["aa_without_glytype"] = True
        elif aa_three in aa_format_dict["glytype"]:
            glycosylation_type = aa_format_dict["glytype"][aa_three]
        else:
            glycosylation_type = ""
            flag_dict["aa_without_glytype"] = True
    elif saccharide.find("|") != -1:
        saccharide_list = saccharide.split("|")
        tmp_list = []
        for sacc in saccharide_list:
            if sacc in glytoucanac2glycosylationtype:
                if "n-linked" in glytoucanac2glycosylationtype[sacc]:
                    tmp_list.append("n-linked")
                if "o-linked" in glytoucanac2glycosylationtype[sacc]:
                    tmp_list.append("o-linked")
                if "c-linked" in glytoucanac2glycosylationtype[sacc]:
                    tmp_list.append("c-linked")
        glycosylation_type = ";".join(tmp_list)
        if glycosylation_type == "" or glycosylation_type == "n-linked;o-linked":
            flag_dict["gtc2glytype_empty_or_n-linked|o-linked"] = True
    elif saccharide in ["G57321FI"]:
        glycosylation_type = "o-linked"
    elif saccharide not in glytoucanac2glycosylationtype:
        glycosylation_type = ""
        if "glycan_not_in_masterlist" not in flag_dict:
            flag_dict["glycan_without_glytype"] = True
    elif glytoucanac2glycosylationtype[saccharide] == []:
        glycosylation_type = ""
    elif amino_acid.find("|") == -1: #resume normal function if we have reported glycosylation information 
        tmp_list = []
        if saccharide in glytoucanac2glycosylationtype:
            if "n-linked" in glytoucanac2glycosylationtype[saccharide]:
                tmp_list.append("n-linked")
            if "o-linked" in glytoucanac2glycosylationtype[saccharide]:
                tmp_list.append("o-linked")
            if "c-linked" in glytoucanac2glycosylationtype[saccharide]:
                tmp_list.append("c-linked")
            if aa_three in aa_format_dict["glytype"]:
                if tmp_list == ["n-linked"] and aa_format_dict["glytype"][aa_three] != "n-linked":
                    flag_dict["n_glycan_aa_mismatch"] = True
                if tmp_list == ["o-linked"] and aa_format_dict["glytype"][aa_three] != "o-linked":
                    flag_dict["o_glycan_aa_mismatch"] = True
        glycosylation_type = ";".join(tmp_list)
        if glycosylation_type == "n-linked;o-linked":
            flag_dict["gtc2glytype_n-linked|o-linked"] = True
    elif amino_acid.find("|") != -1:
        aa_list = amino_acid.split("|")
        if aa_list[0] in aa_format_dict["glytype"]:
            glycosylation_type = aa_format_dict["glytype"][aa_list[0]]
        else:
            glycosylation_type = ""
            flag_dict["aa_without_glytype"] = True

    #remove some flags if glytoucan_type is basecomposition or composition
    skip_list = ["glycan_without_glytype", "gtc2glytype_n-linked|o-linked"]
    if saccharide in gtc_type_dict:
        if gtc_type_dict[saccharide] in ["composition", "basecomposition", "basecomposition;composition"]:
            for flag in skip_list:
                if flag in flag_dict:
                    flag_dict.pop(flag)


    return glycosylation_type


def qc_glyco_fuzzy_position(start_pos, end_pos, peptide, canon, flag_dict):
    
    if start_pos < 0 or start_pos > len(seq_hash[canon]):
        flag_dict["invalid_peptide_start_pos"] = True
    elif end_pos < 0 or end_pos > len(seq_hash[canon]) or end_pos < start_pos:
        flag_dict["invalid_peptide_end_pos"] = True
    else:
        peptide_in_canon = seq_hash[canon][start_pos-1:end_pos]
        if peptide.upper() != peptide_in_canon.upper():
            flag_dict["peptide_mismatch"] = True

    return


def qc_glyco_position(position, amino_acid, canon, flag_dict):

    canon_seq = seq_hash[canon]
    aa_three = amino_acid[0].upper() +  amino_acid[1:].lower() if amino_acid != "" else ""
    if position.find("|") != -1:
        pos_list = position.split("|")
        aa_list = amino_acid.split("|")
        for j in xrange(0, len(pos_list)):
            p = pos_list[j]
            if p.isdigit() == False:
                flag_dict["bad_pos"] = True
            if int(p) == 0 or int(p) >= len(canon_seq):
                flag_dict["bad_pos"] = True
            if j < len(aa_list):
                a = aa_list[j].strip()
                if a != "":
                    if a not in aa_format_dict["one"]:
                        flag_dict["invalid_aa"] = True
                    if aa_format_dict["one"][a] != canon_seq[int(p)-1]:
                        flag_dict["aa_mismatch"] = True
    elif aa_three in aa_format_dict["one"]:
        aa_pos = int(position) if position != "" and position.isdigit() else 0
        if aa_pos == 0 or aa_pos >= len(canon_seq):
            flag_dict["bad_pos"] = True
        elif aa_format_dict["one"][aa_three] != canon_seq[aa_pos-1]:
            flag_dict["aa_mismatch"] = True
    else: #### amino acid is not n or o linked
        flag_dict["invalid_aa"] = True
    return  


def run_global_qc(canon, saccharide, position, amino_acid):
        

    flag_dict = {}
    if saccharide != "" and saccharide not in glycan_dict:
        flag_dict["glycan_not_in_masterlist"] = True

    canon_seq = seq_hash[canon] if canon in seq_hash else ""
    aa_three = amino_acid[0].upper() +  amino_acid[1:].lower() if amino_acid != "" else ""
    if canon_seq != "" and position.find("|") != -1:
        pos_list = position.split("|")
        aa_list = amino_acid.split("|")
        for j in xrange(0, len(pos_list)):
            p = pos_list[j]
            if p.isdigit() == False:
                flag_dict["pos_non_numeric"] = True
            if int(p) == 0 or int(p) > len(canon_seq):
                flag_dict["pos_out_of_seq_range"] = True
            if j < len(aa_list):
                a = aa_list[j].strip()
                if a != "":
                    if a not in aa_format_dict["one"]:
                        flag_dict["invalid_aa"] = True
                    if aa_format_dict["one"][a] != canon_seq[int(p)-1]:
                        flag_dict["aa_mismatch"] = True
    elif aa_three in aa_format_dict["one"]:
        if canon_seq != "" and position.strip() != "":
            aa_pos = int(position) if position.isdigit() else 0
            if aa_pos == 0 or aa_pos > len(canon_seq):
                flag_dict["pos_out_of_seq_range"] = True
            elif aa_format_dict["one"][aa_three] != canon_seq[aa_pos-1]:
                flag_dict["aa_mismatch"] = True
    elif aa_three != "":
        flag_dict["invalid_aa"] = True
   
    glyco_type = qc_glyco_type(amino_acid,saccharide,flag_dict)


    return flag_dict




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
                if xref_id.isdigit() == False:
                    continue
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
   
    longname2shortname = {}
    tax_name_list = []
    for tax_id in species_obj:
        if species_obj[tax_id]["is_reference"] == "yes":
            tax_name = species_obj[tax_id]["long_name"].lower()
            if tax_name not in tax_name_list:
                tax_name_list.append(tax_name)


    canon2genename, genename2canon = load_canon2genename(species)

    out_df = {"fields":[], "data":[]}
    newrow = ["uniprotkb_canonical_ac", "phosphorylation_site_uniprotkb", "amino_acid",
        "kinase_uniprotkb_canonical_ac","kinase_gene_name", 
        "xref_key", "xref_id", "src_xref_key", "src_xref_id"
    ]
    out_df["fields"] = newrow

    FL = open("logs/%s_proteoform_phosphorylation_sites_iptmnet.log" % (species), "w")

    seen_row = {}
    data_frame = {}
    in_file = "downloads/iptmnet/current/ptm.txt"
    log_non_csv_file_use(in_file)


    with open(in_file, "r") as FR:
        for line in FR:
            row = line.split("\t")
            flag_list = []
            if row[0] not in ["PHOSPHORYLATION"]:
                flag_list.append("no-phosphorylation-record")
            tax_name = row[4].split("(")[0].strip().lower()
            if tax_name.find("saccharomyces cerevisiae") != -1 :
                tax_name = "saccharomyces cerevisiae s288c"
            if tax_name.find("drosophila melanogaster") != -1:
                tax_name = "drosophila melanogaster"
 
            #cond_list = []
            #for t_name in tax_name_list:
            #    cond_list.append(tax_name.find(t_name) != -1)
            #if True not in cond_list:
            #    flag_list.append("not-species-of-interest")
            if tax_name != species_obj[species]["long_name"].lower():
                continue
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
            kinase_ac, kinase_canon = row[6].strip(), ""
            if kinase_ac != "" and row[6] in ac2canon:
                kinase_canon = ac2canon[kinase_ac]
            kinase_gene_name = ""
            if kinase_canon in canon2genename:
                kinase_gene_name = canon2genename[kinase_canon]
            xref_key = "protein_xref_" + row[1].lower()
            if  xref_key == "protein_xref_unip":
                xref_key = "protein_xref_uniprotkb"
            xref_id = uniprotkb_ac.split("-")[0]
            #src_xref_key, src_xref_id = "protein_xref_iptmnet", uniprotkb_ac.split("-")[0]
            src_xref_key, src_xref_id = "protein_xref_iptmnet", uniprotkb_ac
            newrow = [canon, str(aa_pos), aa_three, kinase_canon,kinase_gene_name, 
                xref_key,xref_id, src_xref_key, src_xref_id]
            if flag_list != []:
                newrow.append(";".join(flag_list))
                FL.write("\"%s\"\n"  % ("\",\"".join(newrow)))
                continue

            #For original xref_key that comes with the data
            row_str = json.dumps(newrow)
            if row_str not in seen_row:
                out_df["data"].append(newrow)
                seen_row[row_str] = True

            #For xref_key == protein_xref_iptmnet
            xref_id = uniprotkb_ac.split("-")[0]
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
                    if xref_id.isdigit() == False:
                        continue
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
        elif canon.strip() == "":
            flag_list.append("uknown_canon")
        elif canon not in seq_hash:
            flag_list.append("canon_without_seq")
        elif int(aa_pos) == 0 or int(aa_pos) > len(seq_hash[canon]):
            flag_list.append("aa_pos_out_of_range")
        elif aa_format_dict["one"][aa_three] != seq_hash[canon][int(aa_pos)-1]:
            flag_list.append("aa_mismatch")


    return flag_list



def extract_glycosylation_sites_uniprotkb_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()

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
    newrow = required_output_fields 
    newrow += ["uniprotkb_glycosylation_annotation_comment","data_source","eco_id",
        "uniprotkb_ftid", "carb_name", "glycosylation_subtype"]
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
                    out_df["data"].append(newrow)
                    seen_row[newrow_str] = True


    FL.close()

    
    return out_df




def extract_glycation_sites_uniprotkb_ds(species):
    
    data_grid = {"glycosylation":[]}
    sparqlutil.load_glycosylation_sites(data_grid, species)
    out_df = {"fields":[], "data":[]}
    newrow = required_output_fields
    newrow += ["uniprotkb_glycation_annotation_comment","data_source","eco_id","carb_name"]
    for i in xrange(0, len(newrow)):
        if newrow[i].find("glycosylation") != -1:
            newrow[i] = newrow[i].replace("glycosylation", "glycation")
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
                    seen_row[newrow_str] = True


    FL.close()

    return out_df






def extract_glycosylation_sites_pdb_ds(species):

    tissueid2name = load_tissueid2name_map()
    clid2name = load_clid2name_map()


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
    newrow = required_output_fields 
    newrow += ["pdb_id","amino_acid_chain","glycosylation_site_pdb","carb_name"]
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

    FL.close()

    return out_df

def extract_ml_ready_ds(species, dataset):

    cmd = "cp compiled/%s_proteoform_%s.csv unreviewed/" % (species, dataset)
    x = commands.getoutput(cmd)

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


    start_time = time.time()

    global config_obj
    global sparql
    global graph_uri
    global prefixes
    global data_grid
    global species_obj
    global required_output_fields 
    global path_obj
    global seq_hash
    global aa_format_dict
    global ac2canon
    global is_canon
    global ac2canon_strict
    global citation_xref_keys
    global doid2name
    global uberonid2name
    global aa_format_dict
    global glytoucanac2glycosylationtype
    global glycan_dict
    global uckb2glytoucan
    global ds2bco
    global gtc_type_dict
    global nsequon_dict


    species = options.species
    dataset = options.dataset



    citation_xref_keys = ["glycan_xref_pubmed", "protein_xref_pubmed", "protein_xref_doi"]
    config_obj = json.loads(open("conf/config.json", "r").read())
    path_obj = config_obj["pathinfo"]
    ds2bco = json.loads(open("generated/misc/ds2bco.json", "r").read())


    species_obj = {}
    in_file = config_obj["pathinfo"]["misc"]+ "/species_info.csv"
    libgly.load_species_info(species_obj, in_file)
   
    nsequon_dict = load_nsequon_dict()
 
    in_file = path_obj["unreviewed"] + "%s_protein_canonicalsequences.fasta" % (species)
    seq_hash = load_fasta_sequences(in_file)

    #Loading global variables
    ac2canon, is_canon = load_ac2canon(species)
    ac2canon_strict = load_ac2canon_strict(species)
    aa_format_dict = load_aa_format_dict()
    doid2name, uberonid2name = load_disease_info()
    glytoucanac2glycosylationtype = load_glycosylation_type_two()
    glycan_dict = load_glycan_dict()
    uckb2glytoucan = load_uckb2glytoucan()

    gtc_type_dict = {}
    load_glytoucan_type_dict(gtc_type_dict)
    
     



    required_output_fields = ["uniprotkb_canonical_ac","glycosylation_site_uniprotkb",
        "amino_acid","saccharide","glycosylation_type","xref_key",
        "xref_id","src_xref_key","src_xref_id"]


    #start logging file usage
    libgly.log_file_usage("", dataset, "write")


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
    elif dataset == "glycosylation_sites_literature_mining_manually_verified":
        out_df = extract_glycosylation_sites_literature_mining_manually_verified_ds(species)
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
    elif dataset == "glycosylation_sites_oglcnac_mcw":
        out_df = extract_glycosylation_sites_oglcnac_mcw_ds(species)
    elif dataset == "glycosylation_sites_oglcnac_atlas":
        out_df = extract_glycosylation_sites_oglcnac_atlas_ds(species)
    elif dataset in ["glycosylation_sites_o_gluc", "glycosylation_sites_o_gluc_predicted"]:
        out_df = extract_glycosylation_sites_o_gluc_ds(species, dataset)
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
    elif dataset == "glycosylation_sites_unicarbkb_glycomics_study":
        out_df = extract_glycosylation_sites_unicarbkb_glycomics_study_ds(species)
        add_sequon_info(species, out_df)
    elif dataset == "glycosylation_sites_pdc_ccrcc":
        out_df = extract_glycosylation_sites_pdc_ccrcc_ds(species)
    elif dataset == "glycosylation_sites_diabetes_glycomic":
        out_df = extract_glycosylation_sites_diabetes_glycomic_ds(species)
    elif dataset == "glycosylation_sites_embl":
        out_df = extract_glycosylation_sites_embl_ds(species)
    elif dataset in ["ml_ready_pdc_ccrcc", "ml_ready_diabetes_glycomic"]:
        extract_ml_ready_ds(species, dataset)
        exit()
    elif dataset == "glycosylation_sites_carbbank":
        out_df = extract_glycosylation_sites_carbbank_ds(species)
    elif dataset == "glycosylation_sites_platelet":
        out_df = extract_glycosylation_sites_platelet_ds(species)
    elif dataset == "glycosylation_sites_predicted_isoglyp":
        out_df = extract_glycosylation_sites_predicted_isoglyp_ds(species)
    elif dataset == "glycosylation_sites_c_man":
        out_df = extract_glycosylation_sites_c_man_ds(species)
    elif dataset == "glycosylation_sites_glycosmos":
        out_df = extract_glycosylation_sites_glycosmos_ds(species)
    elif dataset == "glycosylation_sites_tablemaker":
        out_df = extract_glycosylation_sites_tablemaker_ds(species)




    append_flag = False
    if dataset.find("glycosylation_sites_") != -1 and dataset.find("citations_glycosylation_sites_") == -1:
        if "start_pos" not in out_df["fields"]:
            append_flag = True
   

    FL = open("logs/%s_proteoform_%s.global.log" % (species, dataset), "w") 
    if append_flag:
        out_df["fields"] += ["start_pos", "end_pos", "start_aa","end_aa", "site_seq"]
    FL.write("\"%s\"\n" % ("\",\"".join(out_df["fields"]  + ["flag_list"])))
    print "\"%s\"" % ("\",\"".join(out_df["fields"]))
    seen_row = {}
    for row in out_df["data"]:
        if append_flag:
            pos = row[out_df["fields"].index("glycosylation_site_uniprotkb")]
            aa_3 = row[out_df["fields"].index("amino_acid")] if "amino_acid" in out_df["fields"] else ""
            aa_1 = aa_format_dict["one"][aa_3] if aa_3 in aa_format_dict["one"] else ""
            row += [pos, pos, aa_3, aa_3, aa_1]
        if dataset.find("glycosylation_sites_") != -1 and dataset.find("citations_") == -1: 
            aa3 = row[out_df["fields"].index("amino_acid")].strip()
            if len(aa3) > 1:
                row[out_df["fields"].index("amino_acid")] = aa3[0].upper() + aa3[1:].lower()
      
        amino_acid, start_aa, start_pos, end_pos = "", "", 0, 0
        if "start_pos" in out_df["fields"]:
            amino_acid = row[out_df["fields"].index("amino_acid")].strip()
            start_pos = row[out_df["fields"].index("start_pos")].strip()
            end_pos = row[out_df["fields"].index("end_pos")].strip()
            start_aa = row[out_df["fields"].index("start_aa")].strip()

        if start_pos == end_pos and amino_acid == "" and start_pos != 0:
            row[out_df["fields"].index("amino_acid")] = start_aa
       
        canon = row[out_df["fields"].index("uniprotkb_canonical_ac")].strip()
        if canon != "" and start_pos == "" and end_pos == "":
            row[out_df["fields"].index("start_pos")] = ""
            row[out_df["fields"].index("end_pos")] = ""
            row[out_df["fields"].index("start_aa")] = ""
            row[out_df["fields"].index("end_aa")] = ""

        #force pipe separated to empty value
        if "glycosylation_site_uniprotkb" in out_df["fields"]:
            pos = row[out_df["fields"].index("glycosylation_site_uniprotkb")]
            if row[out_df["fields"].index("glycosylation_site_uniprotkb")].find("|") != -1:
                s = row[out_df["fields"].index("glycosylation_site_uniprotkb")].split("|")[0]
                e = row[out_df["fields"].index("glycosylation_site_uniprotkb")].split("|")[1]
                row[out_df["fields"].index("glycosylation_site_uniprotkb")] = ""
                start_aa1, end_aa1 = seq_hash[canon][int(s)-1], seq_hash[canon][int(e)-1]
                start_aa3, end_aa3 = "", ""
                if start_aa1 in aa_format_dict["three"]:
                    start_aa3 = aa_format_dict["three"][start_aa1]
                if end_aa1 in aa_format_dict["three"]:
                    end_aa3 = aa_format_dict["three"][end_aa1]
                row[out_df["fields"].index("start_pos")] = s
                row[out_df["fields"].index("end_pos")] = e
                row[out_df["fields"].index("start_aa")] = start_aa3
                row[out_df["fields"].index("end_aa")] = end_aa3
                row[out_df["fields"].index("site_seq")] = seq_hash[canon][int(s)-1:int(e)]

        flag_dict = {}
        if dataset.find("citations_") == -1:
            canon = row[out_df["fields"].index("uniprotkb_canonical_ac")]
            saccharide = row[out_df["fields"].index("saccharide")] if "saccharide" in out_df["fields"] else ""
            aa_three = row[out_df["fields"].index("amino_acid")]
            ll = ["glycosylation_site_uniprotkb","phosphorylation_site_uniprotkb","glyation_site_uniprotkb"]
            pos = ""
            for f in ll:
                if f in out_df["fields"]:
                    pos = row[out_df["fields"].index(f)]
            flag_dict = run_global_qc(canon, saccharide, pos, aa_three)
        flag_list = list(flag_dict.keys())
        if flag_list == []: 
            row_str = json.dumps(row)
            if row_str not in seen_row:
                print "\"%s\"" % ("\",\"".join(row))
            seen_row[row_str] = True
        else:
            FL.write("\"%s\"\n" % ("\",\"".join(row + [";".join(flag_list)])))


    FL.close() 


    
    mol = "proteoform"
    pid = os.getpid()
    src_file = "usage/file_usage.%s.log" % (pid)
    dst_file = "usage/%s_%s_%s.fu.log" % (species, mol, dataset)
    
    cmd = "cat %s |sort -u > %s" % (src_file, dst_file)
    x = commands.getoutput(cmd)
    
    cmd = "rm -f %s" % (src_file)
    x = commands.getoutput(cmd)







if __name__ == '__main__':
        main()

