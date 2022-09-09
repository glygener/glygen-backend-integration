#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict

import libgly




__version__="1.0"
__status__ = "Dev"


def get_counter(counter, prefix):

    if prefix not in counter:
        counter[prefix] = 1
    else:
        counter[prefix] += 1

    zeros = "00000000"
    return prefix.upper() + zeros[0:8-len(str(counter[prefix]))] + str(counter[prefix])



def print_triple(triple, seen):

    if triple not in seen["trpl"]:
        print  triple
        seen["trpl"][triple] = True





def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-v","--ver",action="store",dest="ver",help="data version")


        (options,args) = parser.parse_args()
        for file in ([options.ver]):
            if not (file):
                parser.print_help()
                sys.exit(0)


        config_json = json.loads(open("../conf/config.json", "r").read())
        path_obj  =  config_json[config_json["server"]]["pathinfo"]
        
        ver_dir = "/data/shared/glygen/releases/data/v-%s/" % (options.ver)

        glycan_obj_list = {}
        for in_file in glob.glob(ver_dir + "/jsondb/glycandb/*.json"):
            glytoucan_ac = in_file.split("/")[-1].split(".")[0]
            obj = json.loads(open(in_file, "r").read())
            glycan_obj_list[glytoucan_ac] = obj


        seen = {"trpl":{}, "res":{}, "rxn":{}}
        counter = {}
        ns_map =  config_json["nsmap"]
        uri_map = config_json["urimap"]

        species_obj = {}
        in_file = ver_dir + "/misc/species_info.csv"
        libgly.load_species_info(species_obj, in_file)
        species_list = []
        for k in species_obj:
            obj = species_obj[k]
            if obj["short_name"] not in species_list and obj["is_reference"] == "yes":
                species_list.append(obj["short_name"])



        data_dir = ver_dir + "/reviewed/"
        ac2ec = {}
        for species in species_list:
            in_file = data_dir + "/%s_protein_xref_brenda.csv" % (species)
            if os.path.isfile(in_file) == False:
                continue
            sheet_obj = {}
            libgly.load_sheet_as_dict(sheet_obj, in_file, ",", "uniprotkb_canonical_ac")
            tmp_fl = sheet_obj["fields"]
            for main_id in sheet_obj["data"]:
                for row in sheet_obj["data"][main_id]:
                    ac = main_id.split("-")[0]
                    ac2ec[ac] = row[0]
        
        #generate glycosequence internal ids
        glycanac2seqid = {}
        for glycan_ac in glycan_obj_list:
            for seq_format in ['wurcs', 'glycoct', 'iupac', 'glycam', 'inchi','smiles_isomeric']:
                combo_ac = "%s_%s" % (glycan_ac, seq_format)
                glycanac2seqid[combo_ac] = get_counter(counter, "GLYCOSEQ")



        for glycan_ac in glycan_obj_list:
            obj = glycan_obj_list[glycan_ac]

            glycan_url = uri_map["glycan:Saccharide"] % (glycan_ac)
            triple = "<%s> <%s%s> <%s%s> ." % (glycan_url, ns_map["rdf"], "type", ns_map["glycan"], "Saccharide")
            print_triple(triple, seen)

            image_url = uri_map["glycan:Image"] % (glycan_ac)
            triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"],"has_image", image_url)
            print_triple(triple, seen)

            for key in ['mass']:
                if key in obj:
                    mass_lit = "\"%s\"^^<%s%s>" % (obj[key], ns_map["xsd"], "float")
                    triple = "<%s> <%s%s> %s ." % (glycan_url, ns_map["gly"],"mass", mass_lit)
                    print_triple(triple, seen)

            for seq_format in ['wurcs', 'glycoct', 'iupac', 'glycam', 'inchi','smiles_isomeric']:
                if obj[seq_format] != "":
                    combo_ac = "%s_%s" % (glycan_ac, seq_format)
                    seq_id = glycanac2seqid[combo_ac]
                    seq_url = "%s#%s" % (glycan_url, seq_id)
                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"],"has_glycosequence", seq_url)
                    print_triple(triple, seen)
               
                    triple = "<%s> <%s%s> <%s%s> ." % (seq_url, ns_map["rdf"], "type", ns_map["glycan"], "Glycosequence")
                    print_triple(triple, seen)

                    seq_lit = "\"%s\"" % (obj[seq_format])
                    triple = "<%s> <%s%s> %s ." % (seq_url, ns_map["glycan"], "has_sequence", seq_lit)
                    print_triple(triple, seen)
                
                    format_url = "%s%s%s" % (ns_map["glycan"], "carbohydrate_format_",seq_format)
                    triple = "<%s> <%s%s> <%s> ." % (seq_url, ns_map["glycan"], "in_carbohydrate_format", format_url)
                    print_triple(triple, seen)

            for key in ['species']:
                for o in obj[key]:
                    src_id = get_counter(counter, "SRC")
                    src_url = "%s#%s" % (glycan_url, src_id)
                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "is_from_source", src_url)
                    print_triple(triple, seen)
                    
                    tax_url = "%s%s%s"  % (ns_map["up"], "taxonomy/" , o["taxid"])
                    triple = "<%s> <%s%s> <%s> ." % (src_url, ns_map["glycan"], "has_taxon", tax_url)
                    print_triple(triple, seen)


            for key in ['enzyme']:
                if len(obj[key]) > 0:
                    for o in obj[key]:
                        protein_ac = o["uniprot_canonical_ac"][:-2]
                        uniprot_ac_url = uri_map["up:Protein"] % (protein_ac)
                        ec = ac + "-N/A" 
                        if protein_ac in ac2ec:
                            ec = protein_ac + "-" + ac2ec[protein_ac]
                        if protein_ac not in seen["rxn"]:
                            seen["rxn"][protein_ac] = get_counter(counter, "RXN")
                        rxn_id = seen["rxn"][protein_ac]
                        rxn_url = "%s%s%s"  % (ns_map["gly"], "reaction/" , rxn_id)
                        triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "synthesized_by", rxn_url)
                        print_triple(triple, seen)

                        triple = "<%s> <%s%s> <%s%s> ." % (rxn_url, ns_map["rdf"], "type", ns_map["glycan"], "Glycosyltransferase_Reaction")
                        print_triple(triple, seen)

                        triple = "<%s> <%s%s> <%s> ." % (rxn_url, ns_map["gly"], "has_enzyme_protein", uniprot_ac_url)
                        print_triple(triple, seen)

            for key in ['residues']:
                for o in obj[key]:
                    canon_residue_id = ""
                    if o["id"] not in seen["res"]:
                        canon_residue_id = get_counter(counter, "RES")
                        seen["res"][o["id"]] = canon_residue_id
                    else:
                        canon_residue_id = seen["res"][o["id"]]
                    canon_residue_url = "%s%s%s"  % (ns_map["gly"], "residue/" , canon_residue_id)
                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["gly"], "has_canonical_residue", canon_residue_url)
                    print_triple(triple, seen)
                   
                    if o["parentid"] not in seen["res"]:
                        parent_residue_id = get_counter(counter, "RES")
                        seen["res"][o["parentid"]] = parent_residue_id
                    else:
                        parent_residue_id = seen["res"][o["parentid"]]

                    parent_residue_url = "%s%s%s"  % (ns_map["gly"], "residue/" , parent_residue_id)
                    triple = "<%s> <%s%s> <%s> ." % (canon_residue_url, ns_map["gly"], "has_parent", parent_residue_url)
                    print_triple(triple, seen)



                    lit = "\"%s\"" % (o["name"].lower())
                    triple = "<%s> <%s%s> %s ." % (canon_residue_url, ns_map["gly"], "has_residue_name", lit)
                    print_triple(triple, seen)

                    lit = "\"%s\"" % (o["name"].lower())
                    triple = "<%s> <%s%s> %s ." % (canon_residue_url, ns_map["gly"], "has_residue_id", lit)
                    print_triple(triple, seen)

                    canon = o["attachedby"][4:]
                    if canon not in seen["rxn"]:
                        seen["rxn"][canon] = get_counter(counter, "RXN")
                    rxn_id = seen["rxn"][canon]
                    rxn_url = "%s%s%s"  % (ns_map["gly"], "reaction/" , rxn_id)
                    triple = "<%s> <%s%s> <%s> ." % (canon_residue_url, ns_map["gly"], "attached_by", rxn_url)
                    print_triple(triple, seen)
                    
                    uniprot_ac_url = uri_map["up:Protein"] % (canon[:-2])
                    triple = "<%s> <%s%s> <%s> ." % (rxn_url, ns_map["gly"], "has_enzyme_protein", uniprot_ac_url)
                    print_triple(triple, seen)


            for key in ['crossref']:
                for o in obj[key]:
                    db_url = "%s%s%s" % (ns_map["glycan"], "database_",  o["database"].lower().replace(" ", ""))
                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "glycan_database", db_url)
                    print_triple(triple, seen)

                    label_lit = "\"%s\"" % (o["database"].lower())
                    triple = "<%s> <%s%s> %s ." % (db_url, ns_map["rdfs"], "label", label_lit)
                    print_triple(triple, seen)
        
                    url_lit = "\"%s\"" % (o["url"])
                    triple = "<%s> <%s%s> %s ." % (db_url, ns_map["rdfs"], "recordurl", url_lit)
                    print_triple(triple, seen)
                    
                    urltmpl_lit = "\"%s\"" % ("http://xxx.yy.zzz")
                    triple = "<%s> <%s%s> <%s> ." % (db_url, ns_map["glycan"], "has_url_template", urltmpl_lit)
                    print_triple(triple, seen)

            for key in ["motifs"]:
                for o in obj[key]:
                    motif_ac = o["id"]
                    motif_name = o["name"]
                    motif_url = uri_map["glycan:Glycan_Motif"] % ( motif_ac)
                    triple = "<%s> <%s%s> <%s%s> ." % (motif_url, ns_map["rdf"], "type", 
                            ns_map["glycan"], "Glycan_Motif")
                    print_triple(triple, seen)

                    triple = "<%s> <%s%s> <%s> ." % (glycan_url, ns_map["glycan"], "has_motif", motif_url)
                    print_triple(triple, seen)

                    for seq_format in ['wurcs', 'glycoct', 'iupac', 'glycam', 'inchi','smiles_isomeric']:
                        combo_ac = "%s_%s" % (motif_ac, seq_format)
                        if combo_ac in glycanac2seqid:
                            motif_seq_url = "%s#%s" % (motif_url, glycanac2seqid[combo_ac])
                            triple = "<%s> <%s%s> <%s> ." %(motif_url,ns_map["glycan"],"has_glycosequence", motif_seq_url)
                            print_triple(triple, seen)

                    motif_name_lit = "\"%s\"" % (motif_name)
                    triple = "<%s> <%s%s> %s ." % (motif_url, ns_map["gly"], "has_motif_name", motif_name_lit)
                    print_triple(triple, seen)






if __name__ == '__main__':
        main()




