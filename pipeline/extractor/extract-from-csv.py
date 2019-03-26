import os,sys
import json
import csv

from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment
import commands
import glob



sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"



def extract_ortholog_ds():


    work_book = {}
    xref2canon = {"hgnc":{}, "mgi":{}, "oma":{}}
    canon2xref = {"hgnc":{}, "mgi":{}, "oma":{}}

    xref = "hgnc"
    sheet_obj = {}
    in_file = "unreviewed/human_protein_xref_hgnc.csv"
    libgly.load_sheet(sheet_obj, in_file, ",")
    f_list = sheet_obj["fields"]
    for row in sheet_obj["data"]:
        xref2canon[xref][row[f_list.index("database_id")]] = row[f_list.index("uniprotkb_canonical_ac")]
        canon2xref[xref][row[f_list.index("uniprotkb_canonical_ac")]] = {
            "id":row[f_list.index("database_id")],
            "name":row[f_list.index("database_label")]
        }

    xref = "mgi"
    sheet_obj = {}
    in_file = "unreviewed/mouse_protein_xref_mgi.csv" 
    libgly.load_sheet(sheet_obj, in_file, ",")
    f_list = sheet_obj["fields"]
    for row in sheet_obj["data"]:
        xref2canon[xref][row[f_list.index("database_id")]] = row[f_list.index("uniprotkb_canonical_ac")]
        canon2xref[xref][row[f_list.index("uniprotkb_canonical_ac")]] = {
            "id":row[f_list.index("database_id")],
            "name":row[f_list.index("database_label")]
        }

    xref = "oma"
    for in_file in glob.glob("unreviewed/*_protein_xref_oma.csv"):
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        f_list = sheet_obj["fields"]
        for row in sheet_obj["data"]:
            xref2canon[xref][row[f_list.index("database_id")]] = row[f_list.index("uniprotkb_canonical_ac")]
            canon2xref[xref][row[f_list.index("uniprotkb_canonical_ac")]] = {
                "id":row[f_list.index("database_id")],
                "name":row[f_list.index("database_label")]
            }

    pepid2canon = {}
    sheet_obj = {}
    for in_file in glob.glob("unreviewed/*_protein_transcriptlocus.csv"):
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
    in_file = "downloads/mgi/mouse2human.csv"
    if os.path.isfile(in_file) == True:
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        f_list = sheet_obj["fields"]
        ortho_dict = {}
        homologene2hgnc = {}
        homologene2mgi = {}
        for row in sheet_obj["data"]:
            homologene_id = row[f_list.index("HomoloGene ID")]
            if homologene_id not in ortho_dict:
                ortho_dict[homologene_id] = {}
            if row[f_list.index("NCBI Taxon ID")] == "9606":
                hgnc_id = row[f_list.index("HGNC ID")].split(":")[-1]
                homologene2hgnc[homologene_id] = hgnc_id
                if hgnc_id in xref2canon["hgnc"]:
                    canon = xref2canon["hgnc"][hgnc_id]
                    ortho_dict[homologene_id]["human"] = [canon, row[f_list.index("NCBI Taxon ID")]]
            elif row[f_list.index("NCBI Taxon ID")] == "10090":
                mgi_id = row[f_list.index("Mouse MGI ID")].split(":")[-1]
                homologene2mgi[homologene_id] = mgi_id
                if mgi_id in xref2canon["mgi"]:
                    canon = xref2canon["mgi"][mgi_id]
                    ortho_dict[homologene_id]["mouse"] = [canon, row[f_list.index("NCBI Taxon ID")]]

        for homologene_id in ortho_dict:
            if "human" in ortho_dict[homologene_id] and "mouse" in ortho_dict[homologene_id]:
                obj_list = [
                    {"datasetname":"mgi", "datasetid":homologene2mgi[homologene_id]},
                    {"datasetname":"hgnc", "datasetid":homologene2hgnc[homologene_id]},
                    {"datasetname":"homologene", "datasetid":homologene_id},        
                ]
                for obj in obj_list:
                    row = ortho_dict[homologene_id]["mouse"] + ortho_dict[homologene_id]["human"]
                    row += [obj["datasetname"], obj["datasetid"]]
                    row_string = json.dumps(row)
                    if row_string not in seen_out_row:
                        out_rows.append(row)
                        seen_out_row[row_string] = True



    n1, n2, n3 = 0, 0, 0
    in_file = "downloads/oma/mouse2human.csv"
    if os.path.isfile(in_file) == True:
        peptide_map = {}
        sheet_obj = {}
        libgly.load_sheet(sheet_obj, in_file, ",")
        f_list = sheet_obj["fields"]
        for row in sheet_obj["data"]:
            if row[f_list.index("ortholog_subtype")] == "1:1":
                n1 += 1
                peptide_id = row[f_list.index("ensembl_peptide_id")].split(".")[0]
                peptide_id_ortholog = row[f_list.index("ensembl_peptide_id_ortholog")].split(".")[0]
                peptide_map[peptide_id] = peptide_id_ortholog
                if peptide_id_ortholog in pepid2canon and peptide_id in pepid2canon:
                    n2 += 1
                    pep_one, pep_two = peptide_id, peptide_id_ortholog
                    canon_one, canon_two = pepid2canon[peptide_id], pepid2canon[peptide_id_ortholog]
                    if canon_one in canon2xref["oma"] and canon_two in canon2xref["oma"]:
                        n3 += 1
                        oma_id = canon2xref["oma"][canon_one]["id"]
                        newrow = [canon_one, "10090", canon_two, "9606","oma", oma_id]
                        newrow_string = json.dumps(newrow)
                        if newrow_string not in seen_out_row:
                            out_rows.append(newrow)
                            seen_out_row[newrow_string] = True

    row = [
        "uniprotkb_canonical_ac","tax_id", "uniprotkb_canonical_ac_ortholog", "tax_id_ortholog",
        "evidence_dataset_name", "evidence_dataset_id"
    ]
    print "\"%s\"" % ("\",\"".join(row))
    for row in out_rows:
        print "\"%s\"" % ("\",\"".join(row))



    return





def extract_expression_disease_ds(species, in_file):

    work_book = {}
    sheet_name = "idmapping"
    work_book[sheet_name] = {}
    idmap_file = "unreviewed/%s_protein_idmapping.csv" % (species)
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
    sheet_name = "idmapping"
    work_book[sheet_name] = {}
    idmap_file = "unreviewed/%s_protein_idmapping.csv" % (species)
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






def extract_mutation_ds(species, in_file):



    work_book = {}
    sheet_name = "idmapping"
    work_book[sheet_name] = {}
    idmap_file = "unreviewed/%s_protein_idmapping.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")

    canon_list = []
    for row in work_book[sheet_name]["data"]:
        canon_list.append(row[0])

       
    row = ["uniprotkb_canonical_ac","aa_pos","ref_aa","alt_aa","chromosome_id","chr_pos","ref_nt","alt_nt", "mut_freq","data_source","do_id","do_name"]
    print "\"%s\"" % ("\",\"".join(row))

    f_list = []
    with open(in_file, "r") as FR:
        csv_grid = csv.reader(FR, delimiter=",", quotechar='\"')
        row_count = 0
        for row in csv_grid:
            row_count += 1
            if row_count == 1:
                f_list = row
                continue
            isoform = row[f_list.index("uniprot_ac")]
            if isoform not in canon_list:
                continue
            if int(row[f_list.index("mut_freq")]) < 10:
                continue
            if row[f_list.index("do_id")] in ["3963", "0070003", "3119"]:
                continue
            newrow = [
                isoform,
                row[f_list.index("uniprot_pos")],
                row[f_list.index("ref_aa")],
                row[f_list.index("alt_aa")],
                row[f_list.index("chr_id")],
                row[f_list.index("chr_pos")],
                row[f_list.index("ref_nt")],
                row[f_list.index("alt_nt")],
                row[f_list.index("mut_freq")],
                row[f_list.index("data_src")],
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

    gtf_file = "/data/external/ucsc/grch38/gtf/Homo_sapiens.GRCh38.95.gtf"
    gene_file = "unreviewed/human_protein_xref_hgnc.csv"
    if species == "mouse":
        gtf_file = "/data/external/ucsc/grch38/gtf/Mus_musculus.GRCm38.95.gtf"
        gene_file = "unreviewed/mouse_protein_xref_mgi.csv"


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
                    newrow = [att_dict["gene_id"],row[0],row[3],row[4]]
                    gene_locus[att_dict["gene_name"]] = newrow


    sheet_obj = {}
    libgly.load_sheet(sheet_obj, gene_file, ",")
    newrow = ["uniprotkb_canonical_ac","gene_symbol", "ensembl_gene_id", "chromosome_id", "start_pos", "end_pos"]
    print "\"%s\"" % ("\",\"".join(newrow))
    for row in sheet_obj["data"]:
        canon, gene_name = row[0], row[2]
        if gene_name in gene_locus:
            newrow = [canon, gene_name] + gene_locus[gene_name]
            print "\"%s\"" % ("\",\"".join(newrow))

    
    return



def extract_glycohydrolase_ds(species):

    work_book = {}
    sheet_name = "idmapping"
    work_book[sheet_name] = {}
    idmap_file = "unreviewed/%s_protein_idmapping.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], idmap_file, ",")

    seen_canon = {}
    for row in work_book[sheet_name]["data"]:
        canon = row[0]
        ac = canon.split("-")[0]
        seen_canon[canon] = True


    sheet_name = "glycohydrolase"
    work_book[sheet_name] = {}
    in_file = "downloads/curation/%s_protein_glycohydrolase.csv" % (species)
    libgly.load_sheet(work_book[sheet_name], in_file, ",")

    row = ["uniprotkb_canonical_ac","status","uniprotkb_protein_name","ec_number","cazy","brenda_ec_number","interpro_id","interpro_family_short_name","pfam_ac","pfam_family_name_short"]
    print "\"%s\"" % ("\",\"".join(row))

    f_list = work_book[sheet_name]["fields"]
    for row in work_book[sheet_name]["data"]:
        if row[0] not in seen_canon:
            continue
        #row[f_list.index("aa_pos")],
        newrow = row[0:3]
        j, sep = f_list.index("ec_number"), " "
        ec_list = row[j].strip().split(sep) if row[j].strip() else []

        j, sep = f_list.index("ec_number"), "|"
        brenda_list = row[j].strip().split(sep) if row[j].strip() else []
        
        j, sep = f_list.index("interpro"), "|"
        interpro_idlist = row[j].strip().split(sep) if row[j].strip() else []

        j, sep = f_list.index("interpro_family_short_name"), "|"
        interpro_namelist = row[j].strip().split(sep) if row[j].strip() else []
        
        j, sep = f_list.index("pfam"), "|"
        pfam_idlist = row[j].strip().split(sep) if row[j].strip() else []
                        
        j, sep = f_list.index("pfam_domain"), "|"
        pfam_namelist = row[j].strip().split(sep) if row[j].strip() else []


        print interpro_idlist, interpro_namelist



###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")
    parser.add_option("-d","--dataset",action="store",dest="dataset",help="[idmapping, transcriptlocus]")



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
    global species_info

    species = options.species
    dataset = options.dataset

    species_info = {
        "human":{"taxid":9606, "taxname":"Homo sapiens"}
        ,"mouse":{"taxid":10090, "taxname":"Mus musculus"}    
    }
    config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())

        
    if dataset == "mutation":
        in_file = "downloads/biomuta/biomuta.csv"
        extract_mutation_ds(species, in_file)
    elif dataset == "expression_disease":
        in_file = "downloads/bioxpress/bioxpress_disease.csv"
        extract_expression_disease_ds(species, in_file)
    elif dataset == "expression_normal":
        in_file = "downloads/bioxpress/bioxpress_normal.csv"
        extract_expression_normal_ds(species, in_file)
    elif dataset == "glycohydrolase":
        extract_glycohydrolase_ds(species)
    elif dataset == "ortholog":
        extract_ortholog_ds()
    elif dataset == "genelocus":
        extract_genelocus_ds(species)


















if __name__ == '__main__':
        main()

