import os,sys,glob
import json
import csv
from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--infile",action="store",dest="infile",help="NT input file")
    (options,args) = parser.parse_args()
    for file in ([options.infile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    in_file = options.infile
    config_json = json.loads(open("../conf/config-1.json", "r").read())

    if "human" in in_file:
        out_file = config_json["pathinfo"]["unreviewed"]+"/"+"human_proteoform_glycosylation_sites_uniprot.csv"
    else:
        out_file = config_json["pathinfo"]["unreviewed"]+"/"+"mouse_proteoform_glycosylation_sites_uniprot.csv"

    uniprotAcc_glyAnn = {}
    gly_region = {}
    region_position = {}
    position_site = {}
    gly_comment = {}
    gly_attribution = {}
    att_source1 = {}
    att_source2 = {}
    att_evi = {}
    amino_acids = ["threonine","asparagine","serine","lysine","hydroxylysine","tryptophan","tyrosine","hydroxyproline","cysteine","histidine","isoleucine","valine"]
    with open(in_file, "r") as FR:
        for line in FR:
            row = line.strip().split(' ')
            if row[1] in "<http://purl.uniprot.org/core/annotation>":
                uniprotAcc = row[0].split("/")[-1].replace(">","")
                if uniprotAcc not in uniprotAcc_glyAnn:
                    uniprotAcc_glyAnn[uniprotAcc]=[]
                uniprotAcc_glyAnn[uniprotAcc].append(row[2].split("/")[-1].replace(">",""))

            if row[1] in "<http://purl.uniprot.org/core/range>":
                glyAnn = row[0].split("/")[-1].replace(">","")
                if glyAnn not in gly_region:
                    gly_region[glyAnn]=""
                gly_region[glyAnn] = row[2].split("/")[-1].replace(">","")
            if row[1] in "<http://biohackathon.org/resource/faldo#begin>":
                region = row[0].split("/")[-1].replace(">","")
                if region not in region_position:
                    region_position[region]=""
                region_position[region] = row[2].split("/")[-1].replace(">","")
            if row[1] in "<http://biohackathon.org/resource/faldo#position>":
                position = row[0].split("/")[-1].replace(">","")
                if position not in position_site:
                    position_site[position]=""
                position_site[position]=row[2]

            if row[1] in "<http://www.w3.org/2000/01/rdf-schema#comment>":
                glyAnn = row[0].split("/")[-1].replace(">","")
                if glyAnn not in gly_comment:
                    gly_comment[glyAnn]=[]
                AA = ""
                for item in amino_acids:
                    if item in line:
                        AA=item
                for item in amino_acids:
                    if item in line:
                        AA=item
                if AA == "threonine":   AA="Thr"
                elif AA == "asparagine":    AA="Asn"
                elif AA == "serine":    AA="Ser"
                elif AA == "lysine":    AA="Lys"
                elif AA == "hydroxylysine": AA="Hyl"
                elif AA == "tryptophan":    AA="Trp"
                elif AA == "tyrosine":  AA="Tyr"
                elif AA == "hydroxyproline":    AA="Hyp"
                elif AA == "cysteine":  AA="Cys"
                elif AA == "histidine": AA="His"
                elif AA == "isoleucine":    AA="Ile"
                elif AA == "valine":    AA="Val"
                else: AA=""
                gly_comment[glyAnn] += [row[2].replace('"',""), row[3].replace('(', '').replace(')', ''), AA]

            if row[1] in "<http://purl.uniprot.org/core/evidence>":
                attribution = row[0].split("/")[-1].replace(">","")
                if attribution not in att_evi:
                    att_evi[attribution]=[]
                att_evi[attribution].append(row[2].split("/")[-1].replace(">",""))
            if row[1] in "<http://purl.uniprot.org/core/attribution>":
                glyAnn = row[0].split("/")[-1].replace(">","")
                if glyAnn not in gly_attribution:
                    gly_attribution[glyAnn] = []
                gly_attribution[glyAnn].append(row[2].split("/")[-1].replace(">",""))
            if row[1] in "<http://purl.uniprot.org/core/source>":
                attribution = row[0].split("/")[-1].replace(">","")
                if attribution not in att_source1:
                    att_source1[attribution]=[]
                att_source1[attribution].append(row[2].split("/")[-1].replace(">",""))
                if attribution not in att_source2:
                    att_source2[attribution]=[]
                if "prosite" in row[2]:
                    att_source2[attribution].append("PROSITE")
                if "hamap" in row[2]:
                    att_source2[attribution].append("HAMAP")
                if "pdb" in row[2]:
                    att_source2[attribution].append("PDB")
                if "SIP" in row[2]:
                    att_source2[attribution].append("UniProt Citations")
                elif "http://purl.uniprot.org/citations/" in row[2]:
                    att_source2[attribution].append("PubMed")
                if "http://purl.uniprot.org/uniprot/" in row[2]:
                    att_source2[attribution].append("UniProt")
    gly_site = {}
    for glyAnn in gly_region:
        region = gly_region[glyAnn]
        position = region_position[region]
        gly_site[glyAnn] = position_site[position].split("^^")[0].replace('"',"")

    gly_source1 = {}
    gly_source2 = {}
    gly_evidence = {}
    for glyAnn in gly_attribution:
        attributions = gly_attribution[glyAnn]
        for element in attributions:
            if glyAnn not in gly_source1:
                gly_source1[glyAnn]=[]
                gly_source2[glyAnn]=[]
                gly_evidence[glyAnn]=[]
                if element in att_source1:
                    gly_source1[glyAnn] += att_source1[element]
                if element in att_evi:
                    gly_evidence[glyAnn] += att_evi[element]
                if element in att_source2:
                    gly_source2[glyAnn] += att_source2[element]

    headerList = ["uniprot_canonical_ac", "glycosylation_type", "saccharide", "amino_acid", "glycosylation_site", "source", "evidence", "evidence_tag","evidence_uniprot"]
    with open(out_file,'w') as outfile1:
        writer = csv.writer(outfile1, delimiter=',')
        writer.writerow(headerList)
        for uniprotAcc in uniprotAcc_glyAnn:
            glyAnn = uniprotAcc_glyAnn[uniprotAcc]
            for item in glyAnn:
                row = [uniprotAcc]
                for element in gly_comment[item]:
                    row += [element]
                row += [gly_site[item]]
                if item in gly_source2:
                    row += ["|".join(gly_source2[item])]
                else:
                    row +=[""]
                if item in gly_source1:
                    row += ["|".join(gly_source1[item])]
                else:
                    row += [""]
                if item in gly_evidence:
                    row += ["|".join(gly_evidence[item])]
                else:
                    row += [""]
                if "#" in item:
                    row += [""]
                else:
                    row += [item]
                writer.writerow(row)

if __name__ == '__main__':
    main()
