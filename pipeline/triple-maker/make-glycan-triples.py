#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict


sys.path.append('../../glytools/')
import libgly




__version__="1.0"
__status__ = "Dev"


##########################
def get_counter(counter, prefix):

    if prefix not in counter:
        counter[prefix] = 1
    else:
        counter[prefix] += 1

    zeros = "00000000"
    return prefix.upper() + zeros[0:8-len(str(counter[prefix]))] + str(counter[prefix])



#######################################
def main():

        config_json = json.loads(open("../../conf/config-1.1.json", "r").read())
        path_obj  =  config_json[config_json["server"]]["pathinfo"]

        jsondb_file = path_obj["jsondbpath"] + "/glycandb.json"
        glycan_obj_list = json.loads(open(jsondb_file, "r").read())

        seen_triple = {}
        counter = {}
        ns_map =  config_json["nsmap"]
        uri_map = config_json["urimap"]

        data_dir = "unreviewed/"
        proteinac2ec = {}
        for species in ["human", "mouse"]:
            in_file = data_dir + "/%s_protein_xref_brenda.csv" % (species)
            sheet_obj = {}
            libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "uniprotkb_canonical_ac")
            tmp_fl = sheet_obj["fields"]
            for main_id in sheet_obj["data"]:
                for row in sheet_obj["data"][main_id]:
                    ac = main_id.split("-")[0]
                    proteinac2ec[ac] = row[0]
        
        #generate glycosequence internal ids
        glycanac2seqid = {}
        for glycan_ac in glycan_obj_list:
            glycanac2seqid[glycan_ac] = get_counter(counter, "GLYCOSEQ")



        for glycan_ac in glycan_obj_list:
            obj = glycan_obj_list[glycan_ac]

            glycan_url = uri_map["glycan:Saccharide"] % (glycan_ac)
            triple = "<%s> <%s%s> <%s%s> ." % (glycan_url, ns_map["rdf"], "type", ns_map["glycan"], "Saccharide")
            if triple not in seen_triple:
                print triple
                seen_triple[triple] = True

            image_url = uri_map["glycan:Image"] % (glycan_ac)
            triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"],"has_image", image_url)
            if triple not in seen_triple:
                print triple
                seen_triple[triple] = True

            for key in ['mass']:
                mass_lit = "\"%s\"^^<%s%s>" % (obj[key], ns_map["xsd"], "float")
                triple = "<%s> <%s%s> %s ." % (glycan_url, ns_map["gly"],"mass", mass_lit)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

            for key in ['wurcs', 'glycoct', 'iupac']:
                if obj[key] != "":
                    seq_id = glycanac2seqid[glycan_ac]
                    seq_url = "%s#%s" % (glycan_url, seq_id)
                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"],"has_glycosequence", seq_url)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True
                
                    seq_lit = "\"%s\"^^<%s%s>" % (obj[key], ns_map["xsd"], "string")
                    triple = "<%s> <%s%s> %s ." % (seq_url, ns_map["glycan"], "has_sequence", seq_lit)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True
                
                    format_url = "%s%s%s" % (ns_map["glycan"], "carbohydrate_format_",key)
                    triple = "<%s> <%s%s> <%s> ." % (seq_url, ns_map["glycan"], "in_carbohydrate_format", format_url)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True

            for key in ['species']:
                for o in obj[key]:
                    src_id = get_counter(counter, "SRC")
                    src_url = "%s#%s" % (glycan_url, src_id)
                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "is_from_source", src_url)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True
                    
                    tax_url = "<%s%s%s>"  % (ns_map["up"], "taxonomy/" , o["taxid"])
                    triple = "<%s> <%s%s> <%s> ." % (src_url, ns_map["glycan"], "has_taxon", tax_url)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True


            for key in ['enzyme']:
                if len(obj[key]) > 0:
                    for o in obj[key]:
                        protein_ac = o["uniprot_canonical_ac"][:-2]
                        ec = ac + "-N/A" 
                        if protein_ac in proteinac2ec:
                            ec = protein_ac + "-" + proteinac2ec[protein_ac]
                        rxn_id = get_counter(counter, "RXN")
                        rxn_url = "%s#%s" % (glycan_url, rxn_id)
                        triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "synthesized_by", rxn_url)
                        if triple not in seen_triple:
                           print triple
                           seen_triple[triple] = True

                        enzyme_url = "%s%s%s"  % (ns_map["up"], "enzyme/" ,ec) 
                        triple = "<%s> <%s%s> <%s> ." % (rxn_url, ns_map["glycan"], "has_enzyme", enzyme_url)
                        if triple not in seen_triple:
                            print triple
                            seen_triple[triple] = True


            for key in ['crossref']:
                for o in obj[key]:
                    db_url = "%s%s%s" % (ns_map["glycan"], "database_",  o["database"].lower().replace(" ", ""))
                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "glycan_database", db_url)
                    if triple not in seen_triple:
                        print  triple
                        seen_triple[triple] = True

                    label_lit = "\"%s\"^^<%s%s>" % (o["database"].lower(), ns_map["xsd"], "string")
                    triple = "<%s> <%s%s> %s ." % (db_url, ns_map["rdfs"], "label", label_lit)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True

                    triple = "<%s> <%s%s> <%s> ." % (db_url, ns_map["rdfs"], "recordurl", o["url"])
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True
                    urltmpl_lit = "\"%s\"^^<%s%s>" % ("http://xxx.yy.zzz", ns_map["xsd"], "string")
                    triple = "<%s> <%s%s> <%s> ." % (db_url, ns_map["glycan"], "has_url_template", urltmpl_lit)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True

            for key in ["motifs"]:
                for o in obj[key]:
                    motif_ac = o["id"]
                    motif_name = o["name"]
                    motif_url = uri_map["glycan:Glycan_Motif"] % ( motif_ac)
                    triple = "<%s> <%s%s> <%s%s> ." % (motif_url, ns_map["rdf"], "type", 
                            ns_map["glycan"], "Glycan_Motif")
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True

                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "has_motif", motif_url)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True

                    motif_seq_url = "%s#%s" % (motif_url, glycanac2seqid[motif_ac])
                    triple = "<%s> <%s%s> <%s> ." %(motif_url,ns_map["glycan"],"has_glycosequence", motif_seq_url)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True

                    motif_name_lit = "\"%s\"^^<%s%s>" % (motif_name, ns_map["xsd"], "string")
                    triple = "<%s> <%s%s> %s ." % (motif_url, ns_map["glycan"], "has_motif_name", motif_name_lit)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True






if __name__ == '__main__':
        main()




