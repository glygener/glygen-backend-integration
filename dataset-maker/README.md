Config file:
    The file "conf/config.json" contains relevent paths required by the pipeline


Usage:
The following command generates the dataset file 
human_proteoform_glycosylation_sites_unicarbkb.csv.

    python make-proteoform-dataset.py -s human -d glycosylation_sites_unicarbkb > human_proteoform_glycosylation_sites_unicarbkb.csv



Expected input files:
    
    #unicarbkb input file
    in_file = path_obj["downloads"] + "unicarbkb/human_glycosylation_current.csv" 

    #src of glycosylation_type
    in_file = path_obj["downloads"] + "glytoucan/current/export/classification.tsv"

    #glytoucan master list
    in_file = path_obj["downloads"] + "glytoucan/current/export/glycan_properties.tsv"

    #canonical accession master list
    in_file = path_obj["unreviewed"] + "human_protein_idmapping.csv" 

    #protein sequence master list
    fasta_file = path_obj["unreviewed"] + "human_protein_allsequences.fasta"
    
    #amino acid type dictionary
    in_file = path_obj["misc"] +  "aadict.csv"
    

Log files (output)
    log_file_one = path_obj["logs"] + "human_proteoform_glycosylation_sites_unicarbkb.1.log"
    log_file_two = path_obj["logs"] + "human_proteoform_glycosylation_sites_unicarbkb.2.log"



