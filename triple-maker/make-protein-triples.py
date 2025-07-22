#!/usr/bin/python
import os,sys
import string
from optparse import OptionParser
import csv
import json
import glob
from collections import OrderedDict






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


def add_triple (triple):
    if triple not in seen_triple:
        seen_triple[triple] = True
    return


def get_other_triples(doc, sec_name):

    sec_dict = {
        "glycosylation":{
            "type":"Glycosylation_Site",
            "prd":"glycoylationsite",
            "prefix":"GLYCOSITE"
        },
        "phosphorylation":{
            "type":"Phospho_Site", 
            "prd":"phosphosite", 
            "prefix":"PHOSPHOSITE"
        },
        "glycation":{
            "type":"Glycation_Site", 
            "prd":"glycationsite",
            "prefix":"GLYCATIONSITE"
        },
        "snv":{
            "type":"Snv_Site", 
            "prd":"snvsite",
            "prefix":"SNVSITE"
        }
    }
    site_type = sec_dict[sec_name]["type"]
    site_predicate = sec_dict[sec_name]["prd"]
    site_prefix = sec_dict[sec_name]["prefix"]

    canon = doc["uniprot_canonical_ac"]

    canon_url = uri_map["up:Simple_Sequence"] % (canon)
    if len(doc[sec_name]) > 0:
        for o in doc[sec_name]:
            if "start_pos" in o:
                if o["start_pos"] > 0:
                    loc_id = "SITE_%s_%s_%s" % (canon, o["start_pos"], o["end_pos"])
                    loc_url = uri_map["faldo:Region"] % (loc_id)
                    type_uri = "%s%s" % (ns_map["faldo"], "Region" )
                    prd = "%s%s" % (ns_map["rdf"],"type")
                    triple = "<%s> <%s> <%s> ."%(loc_url,prd,type_uri)
                    add_triple(triple)
                        
                    prd = "%s%s" % (ns_map["gly"],"has_site")
                    triple = "<%s> <%s> <%s> ." % (canon_url,prd,loc_url)
                    add_triple(triple)
                   
                    begin_id = "%s_%s" % (canon, o["start_pos"])
                    begin_url = uri_map["faldo:Begin"] % (begin_id)
                    prd = "%s%s" % (ns_map["faldo"],"begin")
                    triple = "<%s> <%s> <%s> ."%(loc_url,prd,begin_url)                       
                    add_triple(triple)
                        
                    pos_lit = "\"%s\"^^<%s%s>" % (o["start_pos"], ns_map["xsd"], "int")
                    prd = "%s%s" % (ns_map["faldo"],"position")
                    triple = "<%s> <%s> %s ." % (begin_url,prd, pos_lit)
                    add_triple(triple)

                    end_id = "%s_%s" % (canon, o["end_pos"])
                    end_url = uri_map["faldo:End"] % (end_id)
                    prd = "%s%s" % (ns_map["faldo"],"end")
                    triple = "<%s> <%s> <%s> ."%(loc_url,prd,end_url)
                    add_triple(triple)

                    pos_lit = "\"%s\"^^<%s%s>" % (o["end_pos"], ns_map["xsd"], "int")
                    prd = "%s%s" % (ns_map["faldo"],"position")
                    triple = "<%s> <%s> %s ." % (end_url,prd, pos_lit)
                    add_triple(triple)
                 
                    site_id = "%s_%s_%s_%s" % (site_prefix, canon, o["start_pos"], o["end_pos"])
                    ns = "gco" if site_type == "Glycosylation_Site" else "gly"
                    site_url = uri_map[ns+ ":" +site_type] % (site_id)
                    prd =  "%s%s" % (ns_map[ns], site_predicate)
                    triple = "<%s> <%s> <%s> ." % (loc_url, prd, site_url)
                    add_triple(triple)

                    prop = "residue"
                    if prop in o:
                        if o[prop].strip() != "":
                            aa_url = uri_map["glycan:Amino_Acid"] % (o[prop])
                            prd = "%s%s" % (ns_map["rdf"],"type")
                            type_uri = "%s%s" % (ns_map["glycan"], "Amino_Acid")
                            triple = "<%s> <%s> <%s> ."%(aa_url,prd,type_uri)
                            add_triple(triple)

                            prd =  "%s%s" % (ns_map["glycan"], "has_amino_acid")
                            triple = "<%s> <%s> <%s> ." % (site_url, prd, aa_url)
                            add_triple(triple)
                     
                    prop_list = ["kinase_uniprot_canonical_ac", "kinase_gene_name", "comment"]
                    prop_list += ["type", "relation"]
                    prop_list += ["chr_id", "chr_pos", "ref_nt", "alt_nt", "ref_aa", "alt_aa"]
                    prop_list += ["minor_allelic_frequency", "glycoeffect", "keywords"]
                    for prop in prop_list:
                        if prop in o:
                            prd = "%s%s"  % (ns_map["gly"], prop) 
                            val_list = []
                            if type(o[prop]) is list:
                                val_list = o[prop]
                            elif o[prop].strip() != "":
                                val_list = [o[prop].strip()]
                            for v in val_list:
                                triple = "<%s> <%s> \"%s\" ." % (site_url, prd, v)
                                add_triple(triple)
                       
                    prop = "glytoucan_ac"
                    if prop in o:
                        if o[prop].strip() != "":
                            glycan_url = uri_map["glycan:Saccharide"] % (o[prop])
                            prd = "%s%s" % (ns_map["rdf"], "type")
                            type_uri = "%s%s" % (ns_map["glycan"], "Saccharide")
                            triple = "<%s> <%s> <%s> ." % (glycan_url, prd, type_uri)
                            add_triple(triple)
                            prd = "%s%s" % (ns_map["gco"],"has_saccharide")
                            triple = "<%s> <%s> <%s> ." % (site_url, prd, glycan_url)
                            add_triple(triple)


                        
                    prop = "evidence"
                    if prop in o:
                        for oo in o[prop]:
                            if oo["database"] != "PubMed":
                                continue
                            ev_id = "%s_%s" % (oo["database"], oo["id"])
                            ev_url = uri_map["glycan:Evidence"] % (ev_id)
                            prd =  "%s%s" % (ns_map["rdf"],"type")
                            type_uri = "%s%s" % (ns_map["glycan"],"Evidence")
                            triple = "<%s> <%s> <%s> ." % (ev_url,prd,type_uri)
                            add_triple(triple)

                            prd = "%s%s" % (ns_map["glycan"], "has_evidence")
                            triple = "<%s> <%s> <%s> ." % (site_url, prd, ev_url)
                            add_triple(triple)

                            cite_url = uri_map["up:Journal_Citation"] % (ev_id)
                            prd = "%s%s" % (ns_map["rdf"],"type")
                            type_uri = "%s%s" % (ns_map["up"],"Journal_Citation")
                            triple = "<%s> <%s> <%s> ." % (cite_url,prd, type_uri)
                            add_triple(triple)

                            prd = "%s%s" % (ns_map["up"], "citation")
                            triple = "<%s> <%s> <%s> ." % (ev_url, prd, cite_url)
                            add_triple(triple)





    return 




def get_glycosylation_triples(doc):


    canon = doc["uniprot_canonical_ac"]
    protein_len = doc["sequence"]["length"]
    glyset_i = {}
    if len(doc["glycosylation"]) > 0:
        for gly_obj in doc["glycosylation"]:
            gsite_obj = {"uniprot_canonical_ac":canon}
            for k in ["start_pos", "end_pos", "residue", "glytoucan_ac"]:
                if k in gly_obj:
                    gsite_obj[k] = gly_obj[k]
                if "start_pos" not in gsite_obj:
                    gsite_obj["start_pos"], gsite_obj["end_pos"] = 1, protein_len
                for ev_obj in gly_obj["evidence"]:
                    evid = ev_obj["id"]
                    if ev_obj["database"] == "PubMed":
                        if evid not in glyset_i:
                            glyset_i[evid] = []
                        glyset_i[evid].append(gsite_obj)


        siteobj2evid = {}
        for evid in glyset_i:
            site_o_list = glyset_i[evid]
            for site_o in site_o_list:
                siteobj_str = json.dumps([site_o])
                if siteobj_str not in siteobj2evid:
                    siteobj2evid[siteobj_str] = []
                siteobj2evid[siteobj_str].append(evid)
        

        gp_obj_list = []
        gpcounter_dict = {}
        for siteobj_str in siteobj2evid:
            evid_list = siteobj2evid[siteobj_str]
            site_o_list = json.loads(siteobj_str)
            site2glycanlist = {}
            for site_o in site_o_list:
                if "start_pos" not in site_o or "end_pos" not in site_o:
                    continue
                #if site_o["start_pos"] != site_o["end_pos"]:
                #    continue
                if "residue" not in site_o:
                    continue
                site = "%s %s" % (site_o["start_pos"], site_o["residue"])
                if site not in site2glycanlist:
                    site2glycanlist[site] = []
                site2glycanlist[site].append(site_o["glytoucan_ac"])
            if canon not in gpcounter_dict:
                gpcounter_dict[canon] = 0
            gpcounter_dict[canon] += 1
            gpid = "%s-%s" % (canon, gpcounter_dict[canon])
           
            gp_obj = {"gpid":gpid, "sitelist":[], "evidlist":evid_list}
            for site in site2glycanlist:
                gp_obj["sitelist"].append({"site":site, "glycanlist":sorted(set(site2glycanlist[site]))})
            gp_obj_list.append(gp_obj)

        glycanset2id = {}
        for gp_obj in gp_obj_list:
            gp_id = gp_obj["gpid"]
            gp_url = uri_map["gly:Glycoprotein"] % (gp_id)
            triple = "<%s> <%s%s> <%s%s> ." % (gp_url, ns_map["rdf"], "type", ns_map["gly"], "Glycoprotein")
            add_triple(triple)
            canon_url = uri_map["up:Simple_Sequence"] % (canon)
            triple = "<%s> <%s%s> <%s> ." % (gp_url, ns_map["up"], "sequence", canon_url)
            add_triple(triple)
           
            pro_url = "gp_id2pro-of-%s" % (gp_id)
            triple = "<%s> <%s%s> <%s> ." % (gp_url, ns_map["gly"], "has_pro_entry", pro_url)
            add_triple(triple)


            for ev in gp_obj["evidlist"]:
                ev_id = "%s_%s" % ("PubMed", ev) 
                ev_url = uri_map["glycan:Evidence"] % (ev_id)
                prd =  "%s%s" % (ns_map["rdf"],"type")
                type_uri = "%s%s" % (ns_map["glycan"],"Evidence")
                triple = "<%s> <%s> <%s> ." % (ev_url,prd,type_uri)
                add_triple(triple)
        
                prd = "%s%s" % (ns_map["glycan"], "has_evidence")
                triple = "<%s> <%s> <%s> ." % (gp_url, prd, ev_url)
                add_triple(triple)

                cite_url = uri_map["up:Journal_Citation"] % (evid)
                prd = "%s%s" % (ns_map["rdf"],"type")
                type_uri = "%s%s" % (ns_map["up"],"Journal_Citation")
                triple = "<%s> <%s> <%s> ." % (cite_url,prd, type_uri)
                add_triple(triple)
               
                prd = "%s%s" % (ns_map["up"], "citation")
                triple = "<%s> <%s> <%s> ." % (ev_url, prd, cite_url)
                add_triple(triple)
                 
                #journal_name_lit = "\"%s\"^^<%s%s>" % ("xxx", ns_map["xsd"], "string")
                #triple = "<%s> <%s%s> %s ." % (cite_url, ns_map["up"], "name", journal_name_lit)
                #add_triple(triple)

                #journal_title_lit = "\"%s\"^^<%s%s>" % ("xxx", ns_map["xsd"], "string")
                #triple = "<%s> <%s%s> %s ." % (cite_url, ns_map["up"], "title", journal_title_lit)
                #add_triple(triple)

            for o in gp_obj["sitelist"]:
                pos, aa = o["site"].split(" ")[0], o["site"].split(" ")[1]
                gp_site_id = "GLYCOSITE_%s_%s" % (gp_id, pos)
                glysite_url = uri_map["gco:Glycosylation_Site"] % (gp_site_id)
                triple = "<%s> <%s%s> <%s%s> ."%(glysite_url,ns_map["rdf"],"type",ns_map["gco"],"Glycosylation_Site")
                add_triple(triple)

                triple = "<%s> <%s%s> <%s> ."%(gp_url,ns_map["gco"],"glycosylated_at",glysite_url)
                add_triple(triple)
                
                exact_pos_id = "EXACTPOS_%s" % (gp_site_id)
                loc_url = uri_map["faldo:ExactPosition"] % (exact_pos_id)
                triple = "<%s> <%s%s> <%s%s> ."%(loc_url,ns_map["rdf"],"type",ns_map["faldo"],"ExactPosition")
                add_triple(triple)
                
                triple = "<%s> <%s%s> <%s> ."%(glysite_url,ns_map["faldo"],"location",loc_url)
                add_triple(triple)
            
                pos_lit = "\"%s\"^^<%s%s>" % (pos, ns_map["xsd"], "int")
                triple = "<%s> <%s%s> %s ." % (loc_url, ns_map["faldo"],"position", pos_lit)
                add_triple(triple)

                aa_url = uri_map["glycan:Amino_Acid"] % (aa)
                triple = "<%s> <%s%s> <%s%s> ."%(aa_url,ns_map["rdf"],"type",ns_map["glycan"],"Amino_Acid")
                add_triple(triple)

                triple = "<%s> <%s%s> <%s> ." % (loc_url, ns_map["glycan"],"has_amino_acid", aa_url)
                add_triple(triple)


                if "" in o["glycanlist"]:
                    o["glycanlist"].remove("")
               
                if len(o["glycanlist"]) > 1:
                    glycan_list_string = json.dumps(o["glycanlist"])
                    if glycan_list_string not in glycanset2id:
                        glycanset_id = get_counter(counter, "GLYCANSET_" + gp_id)
                        glycanset2id[glycan_list_string] = glycanset_id
                        glycanset_url = uri_map["gly:Saccharide_Set"] % (glycanset_id)

                        triple = "<%s> <%s%s> <%s%s> ." % (glycanset_url, ns_map["rdf"], "type", ns_map["gly"], "Saccharide_Set")
                        add_triple(triple)
                             
                        triple = "<%s> <%s%s> <%s> ." % (glysite_url, ns_map["gly"],"has_saccharide_set", glycanset_url)
                        add_triple(triple)

                        for glycan_id in o["glycanlist"]:
                            glycan_url = uri_map["glycan:Saccharide"] % (glycan_id)
                            triple = "<%s> <%s%s> <%s%s> ." % (glycan_url, ns_map["rdf"], "type", ns_map["glycan"], "Saccharide")
                            add_triple(triple)
                            triple = "<%s> <%s%s> <%s> ." % (glycanset_url, ns_map["gco"],"has_saccharide", glycan_url)
                            add_triple(triple)
                elif len(o["glycanlist"]) == 1:
                    glycan_id = o["glycanlist"][0]
                    glycan_url = uri_map["glycan:Saccharide"] % (glycan_id)
                    triple = "<%s> <%s%s> <%s%s> ." % (glycan_url, ns_map["rdf"], "type", ns_map["glycan"], "Saccharide")
                    add_triple(triple)
                    triple = "<%s> <%s%s> <%s> ." % (glysite_url, ns_map["gco"],"has_saccharide", glycan_url)
                    add_triple(triple)

    return 



#######################################
def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-v","--ver",action="store",dest="ver",help="data version")
        parser.add_option("-b","--batch",action="store",dest="batch",help="1/2/3 ...")



        (options,args) = parser.parse_args()
        for file in ([options.ver, options.batch]):
            if not (file):
                parser.print_help()
                sys.exit(0)

        global ns_map
        global uri_map
        global counter
        global seen_triple
        global evid2counterid
        global canonsite2locid 


        config_json = json.loads(open("../conf/config.json", "r").read())
        path_obj  =  config_json[config_json["server"]]["pathinfo"]
        ver_dir = "/data/shared/glygen/releases/data/v-%s/" % (options.ver)

        counter, evid2counterid, seen_triple, canonsite2locid  = {}, {}, {}, {}
        ns_map =  config_json["nsmap"]
        uri_map = config_json["urimap"]

        batch_size = 10000
        batch_idx = int(options.batch)

        file_list = []
        DEBUG = False
        #DEBUG = True
        if DEBUG:
            ac_list = ["P14210"]
            for ac in ac_list:
                file_list += glob.glob(ver_dir + "/jsondb/proteindb/%s*.json" % (ac)) 
        else:
            file_list += glob.glob(ver_dir + "/jsondb/proteindb/*.json")
            start_idx = batch_size*(batch_idx-1)
            end_idx = start_idx + batch_size 
            end_idx = len(file_list) if end_idx > len(file_list)  else end_idx
            file_list = file_list[start_idx:end_idx]
            

        log_file = "logs/protein.%s.log" % (batch_idx)
        with open(log_file, "w") as FW:
            FW.write("")

        out_file = "generated/sparql/glygen/protein.%s.nt" % (batch_idx)
        #out_file = "junk"
        with open(out_file, "w") as FW:
            FW.write("")


        
        is_writen = {}
        idx = 0
        for in_file in file_list:
            idx += 1
            doc = json.loads(open(in_file, "r").read())
            get_glycosylation_triples(doc)
            get_other_triples(doc, "glycosylation")
            get_other_triples(doc, "phosphorylation")
            get_other_triples(doc, "glycation")
            get_other_triples(doc, "snv")
            with open(out_file, "a", encoding="utf-8") as FA:
                for triple in seen_triple:
                    if triple not in is_writen:
                        FA.write("%s\n" %(triple))
                        is_writen[triple] = True
            if idx%1000 == 0:
                 with open(log_file, "a") as FA:
                    FA.write(" ... processed %s protein objects\n" % (idx))




if __name__ == '__main__':
    main()




