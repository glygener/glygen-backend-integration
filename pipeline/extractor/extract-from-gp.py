import os,sys
import json
import csv


from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 
from Bio import SeqIO




__version__="1.0"
__status__ = "Dev"


def extract_refseq_protein_info_ds(species, tax_id):


    in_file = "downloads/ncbi/refseq/refseq_protein_all_%s.gpff" % (tax_id)
    map_file = "unreviewed/%s_protein_xref_refseq.csv" % (species)

    refseq2canon = {}
    canon2refseq = {}
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


    in_file = "downloads/ncbi/refseq/refseq_protein_all_%s.gpff" % (tax_id)
    map_file = "unreviewed/%s_protein_xref_refseq.csv" % (species)


    refseq2canon = {}
    canon2refseq = {}
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

        
    if dataset == "refseq_protein_info":
        extract_refseq_protein_info_ds(species, species_info[species]["taxid"])
    elif dataset == "function_refseq":
        extract_function_refseq_ds(species, species_info[species]["taxid"])





if __name__ == '__main__':
	main()


