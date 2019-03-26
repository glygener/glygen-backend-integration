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

        jsondb_file = path_obj["jsondbpath"] + "/proteindb.json"
        protein_obj_list = json.loads(open(jsondb_file, "r").read())

        counter = {}
        ns_map =  config_json["nsmap"]
        uri_map = config_json["urimap"]


        glyset_i  = {}
        glyset_ii  = {}
        for protein_ac in protein_obj_list:
            obj = protein_obj_list[protein_ac]
            if len(obj["glycosylation"]) > 0:
                for o in obj["glycosylation"]:
                    gly_obj = {"protein_ac":protein_ac, "position":o["position"],
                            "residue":o["residue"], "glytoucan_ac":o["glytoucan_ac"]}
                    for ev_obj in o["evidence"]:
                        evid = ev_obj["id"]
                        if ev_obj["database"] == "PubMed":
                            if protein_ac not in glyset_i:
                                glyset_i[protein_ac] = {}
                            if evid not in glyset_i[protein_ac]:
                                glyset_i[protein_ac][evid] = []
                            glyset_i[protein_ac][evid].append(gly_obj)
                        #elif ev_obj["database"] == "UniCarbKB":
                        #    if protein_ac not in glyset_ii:
                        #        glyset_ii[protein_ac] = {}
                        #    if evid not in glyset_ii[protein_ac]:
                        #        glyset_ii[protein_ac][evid] = []
                        #    glyset_ii[protein_ac][evid].append(gly_obj)

        

        gpcounter_dict = {}
        gset2evid = {}
        for protein_ac in glyset_i:
            for evid in glyset_i[protein_ac]:
                gset = glyset_i[protein_ac][evid]
                gset_string = "%s|%s" % (protein_ac, json.dumps(gset))
                if gset_string not in gset2evid:
                    gset2evid[gset_string] = []
                gset2evid[gset_string].append(evid)

        gp_obj_list = []
        gpcounter_dict = {}
        for gset_string in gset2evid:
            evid_list = gset2evid[gset_string]
            protein_ac,gset = gset_string.split("|")[0], json.loads(gset_string.split("|")[1])
            site2glycanlist = {}
            for o in gset:
                site = "%s %s" % (o["position"], o["residue"])
                if site not in site2glycanlist:
                    site2glycanlist[site] = []
                site2glycanlist[site].append(o["glytoucan_ac"])
            
            if protein_ac not in gpcounter_dict:
                gpcounter_dict[protein_ac] = 0
            gpcounter_dict[protein_ac] += 1
            gpid = "%s-%s" % (protein_ac, gpcounter_dict[protein_ac])
            
            gp_obj = {"gpid":gpid, "sitelist":[], "evidlist":evid_list}
            for site in site2glycanlist:
                gp_obj["sitelist"].append({"site":site, "glycanlist":sorted(set(site2glycanlist[site]))})
            gp_obj_list.append(gp_obj)


        seen_triple = {}
        glycanset2id = {}


        gpsite2locid = {}
        evid2evidenceid = {}
        for gp_obj in gp_obj_list:
            gp_id = gp_obj["gpid"]
            gp_url = uri_map["gly:Glycoprotein"] % (gp_id)
            triple = "<%s> <%s%s> <%s%s> ." % (gp_url, ns_map["rdf"], "type", ns_map["gly"], "Glycoprotein")
            if triple not in seen_triple:
                print triple
                seen_triple[triple] = True
            canon = "-".join(gp_id.split("-")[0:-1])
            canon_url = uri_map["up:Simple_Sequence"] % (canon)
            triple = "<%s> <%s%s> <%s> ." % (gp_url, ns_map["up"], "sequence", canon_url)
            if triple not in seen_triple:
                print triple
                seen_triple[triple] = True
            
            for evid in gp_obj["evidlist"]:
                if evid not in evid2evidenceid:
                    evid2evidenceid[evid] = get_counter(counter, "GLYCOPROTEIN_EVIDENCE")

                ev_url = uri_map["glycan:Evidence"] % (evid2evidenceid[evid])
                triple = "<%s> <%s%s> <%s%s> ." % (ev_url,ns_map["rdf"],"type",ns_map["glycan"],"Evidence")
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

                triple = "<%s> <%s%s> <%s> ." % (gp_url, ns_map["glycan"], "has_evidence", ev_url)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

                cite_url = uri_map["up:Journal_Citation"] % (evid)
                triple = "<%s> <%s%s> <%s%s> ." % (cite_url,ns_map["rdf"],"type",ns_map["up"],"Journal_Citation")
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True
                
                triple = "<%s> <%s%s> %s ." % (ev_url, ns_map["up"], "citation", cite_url)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True
                 
                journal_name_lit = "\"%s\"^^<%s%s>" % ("xxx", ns_map["xsd"], "string")
                triple = "<%s> <%s%s> %s ." % (cite_url, ns_map["up"], "name", journal_name_lit)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

                journal_title_lit = "\"%s\"^^<%s%s>" % ("xxx", ns_map["xsd"], "string")
                triple = "<%s> <%s%s> %s ." % (cite_url, ns_map["up"], "title", journal_title_lit)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

            for o in gp_obj["sitelist"]:
                glysite_id = get_counter(counter, "SITE")
                glysite_url = uri_map["gco:Glycosylation_Site"] % (glysite_id)
                triple = "<%s> <%s%s> <%s%s> ."%(glysite_url,ns_map["rdf"],"type",ns_map["gco"],"Glycosylation_Site")
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

                triple = "<%s> <%s%s> <%s> ."%(gp_url,ns_map["gco"],"glycosylated_at",glysite_url)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True
                
                pos, aa = o["site"].split(" ")[0], o["site"].split(" ")[1]
                gp_site = "%s_%s" % (gp_id, pos) 
                if gp_site not in gpsite2locid:
                    gpsite2locid[gp_site] = get_counter(counter, "LOC")
                
                loc_url = uri_map["faldo:ExactPosition"] % (gpsite2locid[gp_site])
                triple = "<%s> <%s%s> <%s%s> ."%(loc_url,ns_map["rdf"],"type",ns_map["faldo"],"ExactPosition")
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True
                
                triple = "<%s> <%s%s> <%s> ."%(glysite_url,ns_map["faldo"],"location",loc_url)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True
            
                pos_lit = "\"%s\"^^<%s%s>" % (pos, ns_map["xsd"], "int")
                triple = "<%s> <%s%s> %s ." % (loc_url, ns_map["faldo"],"position", pos_lit)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

                aa_url = uri_map["glycan:Amino_Acid"] % (aa)
                triple = "<%s> <%s%s> <%s%s> ."%(aa_url,ns_map["rdf"],"type",ns_map["glycan"],"Amino_Acid")
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True

                triple = "<%s> <%s%s> <%s> ." % (loc_url, ns_map["glycan"],"has_amino_acid", aa_url)
                if triple not in seen_triple:
                    print triple
                    seen_triple[triple] = True


                
                if "" in o["glycanlist"]:
                    o["glycanlist"].remove("")

                if len(o["glycanlist"]) > 1:
                    glycan_list_string = json.dumps(o["glycanlist"])
                    if glycan_list_string not in glycanset2id:
                        glycanset_id = get_counter(counter, "GLYCANSET")
                        glycanset2id[glycan_list_string] = glycanset_id
                        glycanset_url = uri_map["gly:Saccharide_Set"] % (glycanset_id)

                        triple = "<%s> <%s%s> <%s%s> ." % (glycanset_url, ns_map["rdf"], "type", ns_map["gly"], "Saccharide_Set")
                        if triple not in seen_triple:
                            print triple
                            seen_triple[triple] = True
                             
                        triple = "<%s> <%s%s> <%s> ." % (glysite_url, ns_map["gly"],"has_saccharide_set", glycanset_url)
                        if triple not in seen_triple:
                            print triple
                            seen_triple[triple] = True

                        for glycan_id in o["glycanlist"]:
                            glycan_url = uri_map["glycan:Saccharide"] % (glycan_id)
                            triple = "<%s> <%s%s> <%s%s> ." % (glycan_url, ns_map["rdf"], "type", ns_map["glycan"], "Saccharide")
                            if triple not in seen_triple:
                                print triple
                                seen_triple[triple] = True
                            triple = "<%s> <%s%s> <%s> ." % (glycanset_url, ns_map["gco"],"has_saccharide", glycan_url)
                            if triple not in seen_triple:
                                print triple
                                seen_triple[triple] = True
                elif len(o["glycanlist"]) == 1:
                    glycan_id = o["glycanlist"][0]
                    glycan_url = uri_map["glycan:Saccharide"] % (glycan_id)
                    triple = "<%s> <%s%s> <%s%s> ." % (glycan_url, ns_map["rdf"], "type", ns_map["glycan"], "Saccharide")
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True
                    triple = "<%s> <%s%s> <%s> ." % (glysite_url, ns_map["gco"],"has_saccharide", glycan_url)
                    if triple not in seen_triple:
                        print triple
                        seen_triple[triple] = True



if __name__ == '__main__':
        main()




