import os,sys
import json
import csv
import sparqlutil

from optparse import OptionParser


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC



sys.path.append('../../glytools/')
import libgly



__version__="1.0"
__status__ = "Dev"


def extract_glycosylation_sites_uniprotkb_ds(species):

    data_grid = { "ac2canon":{}, "isoformlist":{},  "glycosylation":[]}
    sparqlutil.load_glycosylation_sites(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)

    row = ["uniprotkb_canonical_ac","glycosylation_type","saccharide","amino_acid","glycosylation_site_uniprotkb","data_source","evidence","eco_id","uniprotkb_ftid"]
    print "\"%s\""  % ("\",\"".join(row))
    for row in data_grid["glycosylation"]:
        if row[0] in data_grid["ac2canon"]:
            row[0] = data_grid["ac2canon"][row[0]]
            print "\"%s\""  % ("\",\"".join(row))
    return
            

def extract_function_uniprot_ds(species):

    data_grid = { "ac2canon":{}, "isoformlist":{},  "function":{}}

    sparqlutil.load_function(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)
    
    row = ["uniprotkb_canonical_ac","database_name", "database_id", "evidence","annotation"]
    print "\"%s\""  % ("\",\"".join(row))

    for ac in data_grid["function"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            ann = data_grid["function"][ac] if ac in data_grid["function"] else ""
            if len(ann) > 0:
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


    



def extract_domap_ds(species):

    data_grid = {
        "mimid2doid":{},
        "doid2icd10cmid":{},
        "doid2icd9cmid":{},
        "doid2keggid":{},
        "doid2meshid":{},
        "doid2umlsid":{},
        "doid2name":{},
        "doid2altname":{},
        "doid2def":{},
    }
    sparqlutil.load_do_mapping(data_grid, species)
    
    
    row = ["do_id", "do_name", "xref_database", "id_in_xref_database"]
    print "\"%s\""  % ("\",\"".join(row))
            
    for do_id in data_grid["doid2name"]:
        for do_name in sorted(set(data_grid["doid2name"][do_id])):
            row = [do_id, do_name]
            for k in ["doid2keggid", "doid2icd10cmid", "doid2icd9cmid", "doid2meshid", "doid2umlsid"]:
                if do_id in data_grid[k]:
                    for db_id in sorted(set(data_grid[k][do_id])):
                        target = k[5:]
                        newrow = row + [target, db_id]
                        print "\"%s\""  % ("\",\"".join(newrow))


    for mim_id in data_grid["mimid2doid"]:
        for do_id in data_grid["mimid2doid"][mim_id]:
            for do_name in sorted(set(data_grid["doid2name"][do_id])):
                newrow = [do_id, do_name, "omim", mim_id]
                print "\"%s\""  % ("\",\"".join(newrow))




def extract_disease_ds(species):

    data_grid = {
        "ac2canon":{}, 
        "isoformlist":{}, 
        "ac2omim":{}, 
        "ac2mondo":{}
    }

    glycophenotype_sheet = {}
    in_file = "downloads/ohsu/human_glycophenotype.csv"
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
                print "\"%s\""  % ("\",\"".join(row))
        if ac in data_grid["ac2omim"]:
            for mim_id in data_grid["ac2omim"][ac]:
                row = [canon, "omim", mim_id, ""]
                print "\"%s\""  % ("\",\"".join(row))


    return




def extract_citations_ds(species):

    data_frame = {}
    in_file = "reviewed/%s_protein_blacklisted_pmids.csv" % (species)
    libgly.load_sheet(data_frame, in_file, ",")
    black_list = []
    for row in data_frame["data"]:
        black_list.append(row[0])
    black_list = sorted(set(black_list))


    data_grid = {"ac2canon":{}, "isoformlist":{}, "citelist":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_citelist(data_grid, species)


    row = ["uniprotkb_canonical_ac","pmid","title","journal_name"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["citelist"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            for obj in data_grid["citelist"][ac]:
                if obj["pmid"] not in black_list:
                    row = [canon,obj["pmid"], obj["journaltitle"], obj["journalname"]]
                    print "\"%s\""  % ("\",\"".join(row))



def extract_canonicalsequences_ds(species):

    data_grid = {"ac2canon":{}, "genename":{}, "recnames":{}, "proteinid":{}, "isoformlist":{}, "isoforminfo":{}, "isoformseq":{}}
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)
    sparqlutil.load_isoformseq(data_grid, species)
    sparqlutil.load_proteinid(data_grid, species)
    sparqlutil.load_recnames(data_grid, species)
    sparqlutil.load_genename(data_grid, species)

    tax_name = config_obj["orgname"][species]
    

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

    tax_name = config_obj["orgname"][species]

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

def extract_information_ds(species):

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

def extract_xref_refseq_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}}
    sparqlutil.load_isoformlist(data_grid, species)

    sheet_obj = {}
    sparqlutil.load_refseq2isoform(sheet_obj, species)

    row = ["uniprotkb_canonical_ac","database_id","database_label"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        if canon in sheet_obj:
            row = [canon, sheet_obj[canon][0], sheet_obj[canon][0]]
            print "\"%s\""  % ("\",\"".join(row))
    return


def extract_xrefs_ds(species, ds_name):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "reactome":{}}
    sparqlutil.load_isoformlist(data_grid, species)

    sheet_obj = {}
    sparqlutil.load_ac2xref(sheet_obj, species, config_obj["xref"][ds_name])

    row = ["uniprotkb_canonical_ac","database_id","database_label"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["ac2canon"]:
        canon = data_grid["ac2canon"][ac]
        if ac in sheet_obj:
            for o in sheet_obj[ac]:
                row = [canon, o["id"], o["label"]]
                print "\"%s\""  % ("\",\"".join(row))






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
        if ac in data_grid[ds_name]:
            obj = data_grid[ds_name][ac]
            full_name = obj["fullname"] if "fullname" in obj else ""
            short_name = obj["shortname"] if "shortname" in obj else ""
            row = [canon, full_name, short_name]
            print "\"%s\""  % ("\",\"".join(row))




    return

def extract_transcriptlocus_ds(species):

    data_grid = {"ac2canon":{}, "isoformlist":{}, "isoforminfo":{}, "locusinfo":{}}
    sparqlutil.load_locusinfo(data_grid, species)
    sparqlutil.load_isoformlist(data_grid, species)
    sparqlutil.load_isoforminfo(data_grid, species)

    row = ["uniprotkb_canonical_ac","uniprotkb_isoform_ac","transcript_id","peptide_id","chromosome_id","start_pos","end_pos"]
    print "\"%s\""  % ("\",\"".join(row))
    for ac in data_grid["isoformlist"]:
        if ac in data_grid["ac2canon"]:
            canon = data_grid["ac2canon"][ac]
            for isoform in data_grid["isoformlist"][ac]:
                if isoform in data_grid["locusinfo"]:
                    o = data_grid["locusinfo"][isoform]
                    row = [canon,isoform,o["trsid"], o["pepid"],
                            str(o["chrid"]),str(o["startpos"]),str(o["endpos"])]
                    print "\"%s\"" % ("\",\"".join(row))




def extract_idmapping_ds(species):

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

    species = options.species
    dataset = options.dataset

    config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())
    

    data_grid = {}
    if dataset == "idmapping":
        extract_idmapping_ds(species)
    elif dataset == "transcriptlocus":
        extract_transcriptlocus_ds(species)
    elif dataset == "recnames":
        extract_recnames_ds(species)
    elif dataset == "altnames":
        extract_altnames_ds(species)
    elif dataset == "ac2pdb":
        extract_ac2pdb_ds(species)
    elif dataset == "information":
        extract_information_ds(species)
    elif dataset == "allsequences":
        extract_allsequences_ds(species)
    elif dataset == "canonicalsequences":
        extract_canonicalsequences_ds(species)
    elif dataset == "citations":
        extract_citations_ds(species)
    elif dataset == "disease":
        extract_disease_ds(species)
    elif dataset == "domap":
        extract_domap_ds(species)
    elif dataset == "function_uniprot":
        extract_function_uniprot_ds(species)
    elif dataset == "glycosylation_sites_uniprotkb":
        extract_glycosylation_sites_uniprotkb_ds(species)
    elif dataset in config_obj["xref"]:
        extract_xrefs_ds(species, dataset)
    elif dataset == "xref_refseq":
        extract_xref_refseq_ds(species)
    elif dataset == "test":
        sparqlutil.load_test(data_grid, species)



















if __name__ == '__main__':
	main()
