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


def make_cell_line_clo_to_cellosaurus_mapping_ds():

    newrow = ["clo_id","cell_line_label","cellosaurus_id"]
    print "\"%s\"" % ("\",\"".join(newrow))

    in_file = "downloads/clo/current/clo.owl"
    clo_id, prev_clo_id = "", ""
    open_flag = False
    with open(in_file, "r") as FR:
        for line in FR:
            if line[0:31] == "    <owl:Class rdf:about=\"&obo;":
                clo_id = line.split("=")[1][6:].split("\"")[0]
                open_flag = True
            if line[0:16] == "    </owl:Class>":
                open_flag = False
            if open_flag and line[0:].find("        <rdfs:label rdf:datatype=") != -1:
                clo_label = line.split(">")[1].split("<")[0]
            if line[0:27] == "        <rdfs:seeAlso>RRID:":
                cellosaurus_id = line.split(">")[1].split("<")[0][5:].strip()
                if clo_id != "":
                    newrow = [clo_id, clo_label, cellosaurus_id]
                    print "\"%s\"" % ("\",\"".join(newrow))
    
    return



def make_taxid2name_ds():

    ref_taxid_dict = {
        "9606":{
            "shortname":"human","longname":"Homo sapiens", "commonname":"Human",
            "ntfile":"uniprot-proteome-homo-sapiens.nt","order":1
        },    
        "10090":{
            "shortname":"mouse","longname":"Mus musculus", "commonname":"Mouse",
            "ntfile":"uniprot-proteome-mus-musculus.nt","order":2
        },     
        "10116":{
            "shortname":"rat","longname":"Rattus norvegicus", "commonname":"Rat",
            
            "ntfile":"uniprot-proteome-rattus-norvegicus.nt","order":3
        },
        "11108":{
            "shortname":"hcv1a","longname":"Hepatitis C virus (genotype 1a, isolate H)",
            "commonname":"HCV-H",    
            "ntfile":"uniprot-proteome-hepatitis-c-virus-1a.nt","order":4
        },
        "11116":{
            "shortname":"hcv1b","longname":"Hepatitis C virus (genotype 1b, isolate Japanese)",
            "commonname":"HCV-Japanese",
            "ntfile":"uniprot-proteome-hepatitis-c-virus-1b.nt","order":5
        }, 
        "694009":{
            "shortname":"sarscov1","longname":"SARS coronavirus (SARS-CoV-1)",
            "commonname":"HCoV-SARS",
            "ntfile":"uniprot-proteome-human-sars-coronavirus.nt","order":6
            }, 
        "2697049":{
            "shortname":"sarscov2","longname":"SARS coronavirus (SARS-CoV-2 or 2019-nCoV)",
            "commonname":"SARS-CoV-2",
            "ntfile":"uniprot-proteome-sars-cov-2.nt","order":7,
        }
    }
   
    
    seen = {}
    data_frame = {}
    in_file = path_obj["downloads"] + "glytoucan/current/export/taxa.tsv"
    libgly.load_sheet(data_frame, in_file, "\t")
    f_list = data_frame["fields"]
    for row in data_frame["data"]:
        tax_id = row[f_list.index("TaxID")]
        seen[tax_id] = True
    
    for tax_id in ref_taxid_dict:
        seen[tax_id] = True

    ncbi_name_dict = {}
    in_file = path_obj["downloads"] +  "/ncbi/taxonomy/current/names.dmp"
    with open(in_file, "r") as FR:
        for line in FR:
            parts = line.strip().split("|")
            tax_id = parts[0].strip()
            if tax_id not in seen:
                continue
            if tax_id not in ncbi_name_dict:
                ncbi_name_dict[tax_id] = {}
            if parts[3].strip() == "scientific name":
                ncbi_name_dict[tax_id]["scientific_name"] = parts[1].strip()
            elif parts[3].strip() == "common name":
                ncbi_name_dict[tax_id]["common_name"] = parts[1].strip()

    newrow = ["tax_id","short_name","long_name","common_name",
            "nt_file","is_reference","sort_order"]
    print "\"%s\"" % ("\",\"".join(newrow))
    tax_id_list = ref_taxid_dict.keys()
    for tax_id in seen:
        if tax_id not in tax_id_list:
            tax_id_list.append(tax_id)


    for tax_id in tax_id_list:
        long_name, common_name = "", "" 
        if tax_id in ncbi_name_dict:
            long_name = ncbi_name_dict[tax_id]["scientific_name"]
            common_name = ncbi_name_dict[tax_id]["common_name"] if "common_name" in ncbi_name_dict[tax_id]  else ""
        if long_name == "" and tax_id in ref_taxid_dict:
            long_name = ref_taxid_dict[tax_id]["longname"]
        if tax_id in ref_taxid_dict:
            if "commonname" in ref_taxid_dict[tax_id]:
                common_name = ref_taxid_dict[tax_id]["commonname"]
        short_name = ref_taxid_dict[tax_id]["shortname"] if tax_id in ref_taxid_dict else ""
        sort_order = ref_taxid_dict[tax_id]["order"] if tax_id in ref_taxid_dict else 1000 
        is_ref = "yes" if tax_id in ref_taxid_dict else "no"
        nt_file = ref_taxid_dict[tax_id]["ntfile"] if tax_id in ref_taxid_dict else ""
        newrow = [tax_id, short_name,long_name,common_name,nt_file,is_ref,str(sort_order)]
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


    data_grid = {}
    if dataset == "taxid2name":
        make_taxid2name_ds()
    elif dataset == "cell_line_clo_to_cellosaurus_mapping":
        make_cell_line_clo_to_cellosaurus_mapping_ds()






if __name__ == '__main__':
        main()


