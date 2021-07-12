import os,sys
import json
import csv

from SPARQLWrapper import SPARQLWrapper, JSON, N3
from rdflib import Graph


sparql = SPARQLWrapper("http://localhost:8890/sparql")
sparql.setReturnFormat(JSON)
    
prefixes = """
    PREFIX faldo: <http://biohackathon.org/resource/faldo#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX gly: <https://sparql.glygen.org/ontology/>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX oboinowl: <http://www.geneontology.org/formats/oboInOwl#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
"""

def load_do_mapping(data_grid):


    graph_uri = "http://sparql.glygen.org#disease"
    qs = prefixes
    qs += " SELECT"
    qs += " ?douri ?doval ?xrefval ?doname ?doaltname ?dodef "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?douri oboinowl:id ?doval . "
    qs += " ?douri rdfs:label ?doname . "
    qs += " optional { ?douri oboinowl:hasDbXref ?xrefval . } "
    qs += " optional { ?douri oboinowl:hasExactSynonym ?doaltname . } "
    qs += " optional { ?douri obo:IAO_0000115 ?dodef . } "
    qs += "}"


    xobj_list_one = [
        {"prefix":"ICD10CM:", "mapname":"doid2icd10cm"},
        {"prefix":"ICD9CM:", "mapname":"doid2icd9cm"},
        {"prefix":"KEGG:", "mapname":"doid2kegg"},
        {"prefix":"MESH:", "mapname":"doid2mesh"},
        {"prefix":"UMLS_CUI:", "mapname":"doid2umls"}
    ]
    xobj_list_two = [
        {"prefix":"OMIM:", "mapname":"mimid2doid"},
    ]
    

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            do_id = result["doval"]["value"].replace("DOID:", "").strip()
            do_name = result["doname"]["value"].strip()
            xref_val, do_altname, do_def = "", "", ""
            if "xrefval" in result:
                xref_val = result["xrefval"]["value"].strip()
            if "doaltname" in result:
                do_altname = result["doaltname"]["value"].strip()
            if "dodef" in result:
                do_def = result["dodef"]["value"].strip()
            for o in xobj_list_one:
                if xref_val[0:len(o["prefix"])] == o["prefix"]:
                    xref_id = xref_val[len(o["prefix"]):]
                    if do_id not in data_grid[o["mapname"]]:
                        data_grid[o["mapname"]][do_id] = []
                    data_grid[o["mapname"]][do_id].append(xref_id)
       
            for o in xobj_list_two:
                if xref_val[0:len(o["prefix"])] == o["prefix"]:
                    xref_id = xref_val[len(o["prefix"]):]
                    if do_id not in data_grid[o["mapname"]]:
                        data_grid[o["mapname"]][xref_id] = []
                    data_grid[o["mapname"]][xref_id].append(do_id)
            if do_id not in data_grid["doid2name"]:
                data_grid["doid2name"][do_id] = []
            data_grid["doid2name"][do_id].append(do_name)

            if do_id not in data_grid["doid2altname"]:
                data_grid["doid2altname"][do_id] = []
            data_grid["doid2altname"][do_id].append(do_altname)
            
            if do_id not in data_grid["doid2def"]:
                data_grid["doid2def"][do_id] = []
            data_grid["doid2def"][do_id].append(do_def)





    return

def load_mutagen_annotation(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?startpos ?endpos ?ecoid ?pmid ?comment ?subs "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdf:type up:Mutagenesis_Annotation . "
    qs += " ?annuri up:range ?rangeuri . "
    qs += "?rangeuri faldo:begin ?beginuri . "
    qs += "?rangeuri faldo:end ?enduri . "
    qs += "?beginuri faldo:position ?startpos . "
    qs += "?enduri faldo:position ?endpos . "
    qs += " optional { ?annuri rdfs:comment ?comment . } "
    qs += " optional { ?annuri up:substitution ?subs . } "
    qs += " optional { "
    qs += " ?annuri gly:attribution ?atturi . "
    qs += " ?atturi up:evidence ?ecoid . "
    qs += " ?atturi up:source ?pmid . "
    qs += "}"
    qs += "}"
    
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            start_pos = int(result["startpos"]["value"])
            end_pos = int(result["endpos"]["value"])
            o = {"startpos":start_pos, "endpos":end_pos}
            for k in ["ecoid", "pmid"]:
                o[k] = result[k]["value"].split("/")[-1] if k in result else ""
            for k in ["subs", "comment"]:
                o[k] = result[k]["value"] if k in result else ""
            if ac not in data_grid["mutagenann"]:
                data_grid["mutagenann"][ac] = []
            data_grid["mutagenann"][ac].append(o)


    return 


def load_mimid2disease_name(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?duri ?mimuri ?dname "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?duri rdfs:seeAlso ?mimuri . "
    qs += " ?duri rdf:type up:Disease . "
    qs += " ?duri skos:prefLabel ?dname . "
    qs += "}"
    

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            d_id = result["duri"]["value"].split("/")[-1]
            mim_id = result["mimuri"]["value"].split("/")[-1]
            d_name = result["dname"]["value"]
            if mim_id not in data_grid["mimid2diseasename"]:
                data_grid["mimid2diseasename"][mim_id] = d_name

    return


def load_signalp_annotation(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?startpos ?endpos ?ecoid ?pmid "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdf:type up:Signal_Peptide_Annotation . "
    qs += " ?annuri up:range ?rangeuri . "
    qs += "?rangeuri faldo:begin ?beginuri . "
    qs += "?rangeuri faldo:end ?enduri . "
    qs += "?beginuri faldo:position ?startpos . "
    qs += "?enduri faldo:position ?endpos . "
    qs += " optional { "
    qs += " ?annuri gly:attribution ?atturi . "
    qs += " ?atturi up:evidence ?ecoid . "
    qs += " ?atturi up:source ?pmid . "
    qs += "}"
    qs += "}"

    comment_hash = {}
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            start_pos = int(result["startpos"]["value"])
            end_pos = int(result["endpos"]["value"])
            o = {"startpos":start_pos, "endpos":end_pos}
            for k in ["ecoid", "pmid"]:
                o[k] = result[k]["value"].split("/")[-1] if k in result else ""
            if ac not in data_grid["signalpann"]:
                data_grid["signalpann"][ac] = []
            data_grid["signalpann"][ac].append(o)            

    return


def load_go_annotation(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?acuri ?gouri ?goterm ?gocaturi ?ecouri ?pmiduri "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?acuri up:classifiedWith ?goannuri . "
    qs += "?acuri rdf:type up:Protein . "
    qs += "?goannuri up:classifiedWith ?gouri . "
    qs += "?gouri rdf:type owl:Class . "
    qs += "?gouri rdfs:label ?goterm . "
    qs += "?gouri gly:goClassification ?gocaturi . "
    qs += " optional { "
    qs += " ?goannuri gly:attribution ?atturi . "
    qs += " ?atturi up:evidence ?ecouri . "
    qs += " ?atturi up:source ?pmiduri . "
    qs += "}"
    qs += "}"


    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit) 
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n 
        for result in results["results"]["bindings"]:
            ac = result["acuri"]["value"].split("/")[-1].strip()
            go_id = result["gouri"]["value"].split("/")[-1].strip()
            go_cat = result["gocaturi"]["value"].split("/")[-1].strip()
            go_term = result["goterm"]["value"]
            ecoid, pmid = "", ""
            if "ecouri" in result:
                ecoid = result["ecouri"]["value"].split("/")[-1]
            if "pmiduri" in result:
                pmid = result["pmiduri"]["value"].split("/")[-1]
            if ac not in data_grid["goann"]:
                data_grid["goann"][ac] = []
            o = {"goid":go_id, "goterm":go_term, "gocat":go_cat, 
                    "ecoid":ecoid,"pmid":pmid}
            data_grid["goann"][ac].append(o)



    return



def load_disease_mim(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?xrefuri ?xrefcomment "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?ac rdfs:seeAlso ?xrefuri . "
    qs += "?xrefuri rdfs:comment ?xrefcomment . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit) 
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n 
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1].strip()
            xref_id = result["xrefuri"]["value"].split("/")[-1].strip()
            xref_type = result["xrefuri"]["value"].split("/")[-2].strip()
            xref_comment = result["xrefcomment"]["value"]
            if xref_type == "mim" and xref_comment == "phenotype":
                if ac not in data_grid["ac2omim"]:
                    data_grid["ac2omim"][ac] = []
                data_grid["ac2omim"][ac].append(xref_id)



    return




def load_citelist(data_grid, species):

    pmid2authlist = {}

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?citeuri ?author "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?citeuri up:author ?author . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            pm_id = result["citeuri"]["value"].split("/")[-1]
            author = result["author"]["value"]
            if pm_id not in pmid2authlist:
                pmid2authlist[pm_id] = []
            pmid2authlist[pm_id].append(author)


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?citeuri ?journalname ?journaltitle ?pubdate "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?ac up:citation ?citeuri . "
    qs += "?citeuri up:name ?journalname . "
    qs += "?citeuri up:title ?journaltitle . "
    qs += "?citeuri up:date ?pubdate . "
    qs += "}"

    seen_pmid = {}

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1].strip()
            pm_id = result["citeuri"]["value"].split("/")[-1]
            combo_id = "%s %s" % (ac, pm_id)
            if combo_id not in seen_pmid:
                seen_pmid[combo_id] = True
                journal_name = result["journalname"]["value"].strip()
                pub_date = result["pubdate"]["value"].strip().split("^")[0]
                journal_title = result["journaltitle"]["value"].strip().replace("\"", "`")
                if ac not in data_grid["citelist"]:
                    data_grid["citelist"][ac] = []

                auth_list = pmid2authlist[pm_id] if pm_id in pmid2authlist else []
                o = {"pmid":pm_id, "journalname":journal_name, "journaltitle":journal_title, 
                        "pubdate":pub_date, "authorlist":auth_list}
                data_grid["citelist"][ac].append(o)
    return

                                     


def load_gene_names(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?preflabel ?altlabel ?orfname "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac rdf:type up:Protein . "
    qs += " ?ac up:encodedBy ?geneuri . "
    qs += " optional { ?geneuri skos:prefLabel ?preflabel . } "
    qs += " optional { ?geneuri skos:altLabel ?altlabel . } "
    qs += " optional { ?geneuri up:orfName ?orfname . } "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1].split("#")[0]
            o = {}
            for k in ["preflabel", "altlabel", "orfname"]:
                o[k] = result[k]["value"].split("/")[-1] if k in result else ""
            if ac not in data_grid["genenames"]:
                data_grid["genenames"][ac] = []
            data_grid["genenames"][ac].append(o)

    return


def load_genename(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?genename "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?ac rdf:type up:Protein . "
    qs += "?ac up:encodedBy ?geneuri . "
    qs += "?geneuri skos:prefLabel ?genename . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            gene_name = result["genename"]["value"].split("/")[-1]
            data_grid["genename"][ac] = gene_name

    return


def load_enzymelist_one(data_grid, species):

    ec2activity = {}
    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ecuri ?activity"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?ecuri up:activity ?activityuri . "
    qs += "?activityuri rdfs:label ?activity . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ec = result["ecuri"]["value"].split("/")[-1]
            activity = result["activity"]["value"]
            ec2activity[ec] = activity

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?acuri ?ecuri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?acuri up:enzyme ?ecuri . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["acuri"]["value"].split("/")[-1].split("#")[0]
            ec = result["ecuri"]["value"].split("/")[-1]
            activity = ec2activity[ec] if ec in ec2activity else ""
            o =  {"ec":ec, "activity":activity}
            if ac not in data_grid["enzymeinfo"]:
                data_grid["enzymeinfo"][ac] = []
            data_grid["enzymeinfo"][ac].append(o)

    return


def load_enzymelist_two(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?acuri ?ec"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?acuri up:recommendedName ?nameuri . "
    qs += "?nameuri up:ecName ?ec . "
    qs += "}"


    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["acuri"]["value"].split("/")[-1]
            ec = result["ec"]["value"]
            activity = ""
            o =  {"ec":ec, "activity":activity}
            if ac not in data_grid["enzymeinfo"]:
                data_grid["enzymeinfo"][ac] = []
            data_grid["enzymeinfo"][ac].append(o)
    return


def load_interaction(data_grid, species):

    ac2intact = {}
    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?puri ?puniprotac"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    #qs += "?acuri up:interaction ?intacturi . "
    qs += "?intacturi up:participant ?puri . "
    qs += "?puri owl:sameAs ?puniprotac . "
    #qs += "?intacturi rdf:type up:Self_Interaction . "
    qs += "}"
   

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            pintactac = result["puri"]["value"].split("/")[-1]
            puniprotac = result["puniprotac"]["value"].split("/")[-1]
            puniprotac_part = puniprotac.split("-")[0]
            if puniprotac_part in data_grid["ac2canon"]:
                if data_grid["ac2canon"][puniprotac_part] == puniprotac:
                    puniprotac = puniprotac_part
            if puniprotac not in ac2intact:
                ac2intact[puniprotac] = []
            if pintactac not in ac2intact[puniprotac]:
                ac2intact[puniprotac].append(pintactac)
   

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?acuri ?puri ?intacturi ?experiments ?puniprotac ?pgenename ?puniprotid ?ptaxid "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?acuri up:interaction ?intacturi . "
    qs += "?intacturi up:experiments ?experiments . "
    qs += "?intacturi up:participant ?puri . "
    qs += "?puri owl:sameAs ?puniprotac . "
    qs += "?puri rdfs:label ?pgenename . "
    qs += "?puri up:mnemonic ?puniprotid . "
    qs += "?puri up:organism ?ptaxid . "
    qs += "}"


    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["acuri"]["value"].split("/")[-1]
            if ac not in ac2intact:
                continue
            if len(ac2intact[ac]) > 1:
                continue
            intactac = ac2intact[ac][0]
            pintactac = result["puri"]["value"].split("/")[-1]
            experiments = result["experiments"]["value"]
            puniprotac = result["puniprotac"]["value"].split("/")[-1]
            puniprotid = result["puniprotid"]["value"]
            pgenename = result["pgenename"]["value"]
            ptaxid = result["ptaxid"]["value"].split("/")[-1]
            o =  {"experiments":experiments, "puniprotac":puniprotac, "puniprotid":puniprotid,
                    "pgenename":pgenename, "ptaxid":ptaxid, "intactac":intactac, "pintactac":pintactac}
            if ac not in data_grid["interaction"]:
                data_grid["interaction"][ac] = []
            data_grid["interaction"][ac].append(o)
    return



def load_isoforminfo(data_grid, species):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?isoform ?status ?iscanon"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?isoform rdf:type up:Simple_Sequence . "
    qs += "?isoform up:reviewed ?status . "
    qs += "?isoform gly:canonical ?iscanon . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            isoform = result["isoform"]["value"].split("/")[-1]
            status = result["status"]["value"].split("/")[-1]
            iscanon = result["iscanon"]["value"].split("/")[-1]
            data_grid["isoforminfo"][isoform] = {"reviewed":status, "iscanon":iscanon}

    return


def load_isoformlist(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT" 
    qs += "?ac ?isoform ?iscanon"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?ac rdf:type up:Protein . "
    qs += "?ac up:sequence ?isoform . "
    qs += "?isoform gly:canonical ?iscanon . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            isoform = result["isoform"]["value"].split("/")[-1]
            iscanon = result["iscanon"]["value"].split("/")[-1]
            if ac not in data_grid["isoformlist"]:
                data_grid["isoformlist"][ac] = []
            data_grid["isoformlist"][ac].append(isoform)
            if iscanon == "1":
                data_grid["ac2canon"][ac] = isoform


    return




def load_locusinfo(data_grid, species):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?isoformuri ?trsid ?pepid ?chrid ?startpos ?endpos ?strand "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?isoformuri gly:ensTranscript ?trsid . "
    qs += "?trsid up:translatedTo ?pepid . "
    qs += "?trsid gly:chromosome ?chrid . "
    qs += "?trsid gly:transcriptRange ?rangeuri . "
    qs += "?trsid gly:reverseStrand ?strand . "
    qs += "?rangeuri faldo:begin ?beginuri . "
    qs += "?rangeuri faldo:end ?enduri . "
    qs += "?beginuri faldo:position ?startpos . "
    qs += "?enduri faldo:position ?endpos . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            isoform = result["isoformuri"]["value"].split("/")[-1]
            if isoform not in data_grid["locusinfo"]:
                data_grid["locusinfo"][isoform] = {
                    "chrid":result["chrid"]["value"].strip()
                    ,"trsid":result["trsid"]["value"].split("/")[-1] 
                    ,"pepid":result["pepid"]["value"].split("/")[-1]
                    ,"strand":"x"
                    ,"startpos":result["startpos"]["value"]
                    ,"endpos":result["endpos"]["value"]
                    ,"strand":result["strand"]["value"]
                }

    return



def load_protein_names(data_grid, species, predicate_type):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?nameuri ?fullname ?shortname ?ecname "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac %s ?nameuri ." % (predicate_type)
    #qs += " ?ac rdf:type up:Protein ."
    qs += " ?nameuri up:fullName ?fullname ."
    qs += " optional { ?nameuri up:shortName ?shortname . } "
    qs += " optional { ?nameuri up:ecName ?ecname .  } "
    qs += "}"

    

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1].split("#")[0]
            o = {}
            for k in ["fullname", "shortname", "ecname"]:
                #o[k] = result[k]["value"].split("/")[-1] if k in result else ""
                o[k] = result[k]["value"] if k in result else ""
            if ac not in data_grid["proteinnames"]:
                data_grid["proteinnames"][ac] = []
            data_grid["proteinnames"][ac].append(o)


    return





def load_glycosylation_sites(data_grid, species):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    #qs += " ?ac ?comment ?startpos ?endpos ?evidenceuri ?citationuri"
    qs += " ?ac ?comment ?annuri ?startpos ?endpos ?evidenceuri ?citationuri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdf:type up:Glycosylation_Annotation . "
    qs += " ?annuri rdfs:comment ?comment . "
    qs += " ?annuri up:range ?rangeuri . "
    qs += " ?rangeuri faldo:begin ?beginuri . "
    qs += " ?rangeuri faldo:end ?enduri . "
    qs += " ?beginuri faldo:position ?startpos . "
    qs += " ?enduri faldo:position ?endpos . "
    qs += " optional { "
    qs += " ?annuri gly:attribution ?atturi . "
    qs += " ?atturi up:source ?citationuri . "
    qs += " ?atturi up:evidence ?evidenceuri . "
    qs += "}"
    qs += "}"

    
    aa_dict = {
            "Alanine":"Ala",
            "Arginine":"Arg",
            "Asparagine":"Asn",
            "Aspartic":"acid",
            "Cysteine":"Cys",
            "Glutamic acid":"Gly",
            "Glutamine":"Gln",
            "Glycine":"Gly",
            "Histidine":"His",
            "Isoleucine":"Ile",
            "Leucine":"Leu",
            "Lysine":"Lys",
            "Methionine":"Met",
            "Phenylalanine":"Phe",
            "Proline":"Pro",
            "Serine":"Ser",
            "Threonine":"Thr",
            "Tryptophan":"Trp",
            "Tyrosine":"Tyr",
            "Valine":"Val"
    }


    source_dict = {
        "citations":"PubMed",
        "pdb":"PDB",
        "uniprot":"UniProtKB"
    }
    ev_dict = {}
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            car_id =  result["annuri"]["value"].split("/")[-1]
            car_id = car_id if car_id[0:4] == "CAR_" else ""
            comment = result["comment"]["value"] if "comment" in result else ""
            start_pos = result["startpos"]["value"]
            end_pos = result["endpos"]["value"]
            car_id = car_id if car_id[0:4] == "CAR_" else ""
            evidence, source_one, source_two = "", "", ""
            if "evidenceuri" in result:
                evidence = result["evidenceuri"]["value"].split("/")[-1]
            if "citationuri" in result:
                source_one = result["citationuri"]["value"].split("/")[-2]
                source_two = result["citationuri"]["value"].split("/")[-1]
            gly_type, saccharide = "", ""
            if comment != "":
                gly_type = comment.strip().split(" ")[0]
                saccharide = comment.strip().split(" ")[1]
                saccharide = saccharide.replace("(", "").replace(")", "")
            source_one = source_dict[source_one] if source_one in source_dict else source_one
            for aa in aa_dict:
                if comment.find(aa.lower()) != -1:
                    aa = aa_dict[aa]
                    row = [ac, gly_type, saccharide, aa, start_pos, comment]
                    combo_id = ",".join(row)
                    if combo_id not in ev_dict:
                        ev_dict[combo_id] = []
                    ev_dict[combo_id].append([source_one, source_two,evidence, car_id])

    for combo_id in ev_dict:
        row_one = combo_id.split(",")
        for row_two in ev_dict[combo_id]:
            data_grid["glycosylation"].append(row_one + row_two)
    return





def load_phosphorylation_sites(data_grid, species):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?comment ?annuri ?startpos ?endpos ?evidenceuri ?citationuri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdfs:comment ?comment . "
    qs += " ?annuri rdf:type up:Modified_Residue_Annotation . "
    qs += " ?annuri gly:attribution ?atturi . "
    qs += " ?atturi up:source ?citationuri . "
    qs += " ?atturi up:evidence ?evidenceuri . "
    qs += " ?annuri up:range ?rangeuri . "
    qs += "?rangeuri faldo:begin ?beginuri . "
    qs += "?rangeuri faldo:end ?enduri . "
    qs += "?beginuri faldo:position ?startpos . "
    qs += "?enduri faldo:position ?endpos . "
    qs += "}"

   
    aa_dict = {
            "Alanine":"Ala",
            "Arginine":"Arg",
            "Asparagine":"Asn",
            "Aspartic":"acid",
            "Cysteine":"Cys",
            "Glutamic":"acid",
            "Glutamine":"Gln",
            "Glycine":"Gly",
            "Histidine":"His",
            "Isoleucine":"Ile",
            "Leucine":"Leu",
            "Lysine":"Lys",
            "Methionine":"Met",
            "Phenylalanine":"Phe",
            "Proline":"Pro",
            "Serine":"Ser",
            "Threonine":"Thr",
            "Tryptophan":"Trp",
            "Tyrosine":"Tyr",
            "Valine":"Val"
    }


    source_dict = {
        "citations":"PubMed",
        "pdb":"PDB",
        "uniprot":"UniProtKB"
    }
    ev_dict = {}
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            comment = result["comment"]["value"]
            evidence = result["evidenceuri"]["value"].split("/")[-1]
            if comment[0:7] != "Phospho":
                continue
        

            gly_type = comment.strip().split(" ")[0]
            #saccharide = comment.strip().split(" ")[1]
            #saccharide = saccharide.replace("(", "").replace(")", "")
            start_pos = result["startpos"]["value"]
            end_pos = result["endpos"]["value"]
            source_one = result["citationuri"]["value"].split("/")[-2]
            source_one = source_dict[source_one] if source_one in source_dict else source_one
            source_two = result["citationuri"]["value"].split("/")[-1]

            for aa in aa_dict:
                if comment.find(aa.lower()) != -1:
                    aa = aa_dict[aa]
                    row = [ac, gly_type, aa, start_pos]
                    combo_id = ",".join(row)
                    if combo_id not in ev_dict:
                        ev_dict[combo_id] = []
                    ev_dict[combo_id].append([source_one, source_two,evidence])


    
    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?comment ?annuri ?startpos ?endpos ?isoform "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdfs:comment ?comment . "
    qs += " ?annuri rdf:type up:Modified_Residue_Annotation . "
    qs += " ?annuri up:range ?rangeuri . "
    qs += "?rangeuri faldo:begin ?beginuri . "
    qs += "?rangeuri faldo:end ?enduri . "
    qs += "?beginuri faldo:position ?startpos . "
    qs += "?enduri faldo:position ?endpos . "
    qs += " optional { ?annuri up:sequence ?isoform . } "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            comment = result["comment"]["value"]
            isoform = result["isoform"]["value"] if "isoform" in result else ""
            if comment[0:7] != "Phospho":
                continue

            #ignore annotation given to isoforms
            if isoform != "":
                continue
            kinase_list = [""]
            if comment.find("; by") != -1:
                for gene_word in comment.split("; by")[1].replace(",", "").split(" "):
                    if gene_word.strip() != "":
                        kinase_list += [gene_word.replace(";", "")]
                for w in ["", "kinase", "form", "vitro", 
                        "and", "in", "by", "autocatalysis", "or"]:
                    if w in kinase_list:
                        kinase_list.remove(w)
            gly_type = comment.strip().split(" ")[0].replace(";", "")
            #saccharide = comment.strip().split(" ")[1]
            #saccharide = saccharide.replace("(", "").replace(")", "")
            start_pos = result["startpos"]["value"]
            end_pos = result["endpos"]["value"]
            for aa in aa_dict:
                if comment.find(aa.lower()) != -1:
                    aa = aa_dict[aa]
                    row_one = [ac, gly_type, aa, start_pos]
                    combo_id = ",".join(row_one)
                    if combo_id in ev_dict:
                        for row_two in ev_dict[combo_id]:
                            row_three = row_one + [comment] + row_two
                            data_grid["phosphorylation"].append(row_three)
                            #for gene_name in kinase_list:
                            #    row_three = row_one + [gene_name] + row_two
                            #    data_grid["phosphorylation"].append(row_three)
                    else:
                        row_three = row_one + [comment, "", "", ""]
                        data_grid["phosphorylation"].append(row_three)
                        #for gene_name in kinase_list:
                        #    row_three = row_one + [gene_name, "", "", ""]
                        #    data_grid["phosphorylation"].append(row_three)
    
    return




def load_pdbinfo(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?pdburi ?method ?resolution "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac rdfs:seeAlso ?pdburi . "
    qs += " ?pdburi up:method ?method . "
    qs += " ?pdburi up:resolution ?resolution . "
    qs += " ?pdburi rdf:type up:Structure_Resource . "
    qs += "}"

    comment_hash = {}
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            pdb_id = result["pdburi"]["value"].split("/")[-1]
            method = result["method"]["value"].split("/")[-1]
            res = float(result["resolution"]["value"])
            if ac not in data_grid["pdbinfo"]:
                data_grid["pdbinfo"][ac] = []
            o = {"pdbid":pdb_id, "method":method, "res":res}
            data_grid["pdbinfo"][ac].append(o)

    return



def load_ptm_annotation(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?startpos ?endpos ?ecoid ?pmid ?comment "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdf:type up:PTM_Annotation . "
    qs += " optional { ?annuri rdfs:comment ?comment . } "
    qs += " optional { "
    qs += " ?annuri gly:attribution ?atturi . "
    qs += " ?atturi up:evidence ?ecoid . "
    qs += " ?atturi up:source ?pmid . "
    qs += "}"
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            o = {}
            for k in ["ecoid", "pmid"]:
                o[k] = result[k]["value"].split("/")[-1] if k in result else ""
            for k in ["comment"]:
                o[k] = result[k]["value"] if k in result else ""
            for k in ["startpos", "endpos"]:
                o[k] = int(result[k]["value"]) if k in result else -1

            if ac not in data_grid["ptmann"]:
                data_grid["ptmann"][ac] = []
            data_grid["ptmann"][ac].append(o)
    
      

    return


def load_site_annotation(data_grid, species, ann_type):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?comment ?annuri ?startpos ?endpos ?subs ?pmid ?ecoid "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdfs:comment ?comment . "
    qs += " ?annuri rdf:type up:%s . " % (ann_type)
    qs += " ?annuri up:range ?rangeuri . "
    qs += "?rangeuri faldo:begin ?beginuri . "
    qs += "?rangeuri faldo:end ?enduri . "
    qs += "?beginuri faldo:position ?startpos . "
    qs += "?enduri faldo:position ?endpos . "
    qs += " optional { ?annuri up:substitution ?subs . } "
    qs += " optional { "
    qs += " ?annuri gly:attribution ?atturi . "
    qs += " ?atturi up:evidence ?ecoid . "
    qs += " ?atturi up:source ?pmid . "
    qs += "}"
 
    
    qs += "}"

    ev_dict = {}
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            comment = result["comment"]["value"]
            start_pos = result["startpos"]["value"]
            end_pos = result["endpos"]["value"]
            o = {"anntype":ann_type, "start":start_pos, "end":end_pos, "ann":comment}
            for k in ["ecoid", "pmid"]:
                o[k] = result[k]["value"].split("/")[-1] if k in result else ""
            for k in ["subs"]:
                o[k] = result[k]["value"] if k in result else ""
            if ac not in data_grid["siteann"]:
                data_grid["siteann"][ac] = []
            data_grid["siteann"][ac].append(o)



    return






def load_function(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?ann"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdfs:comment ?ann . "
    qs += " ?annuri rdf:type up:Function_Annotation . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            ann = result["ann"]["value"]
            if ac not in data_grid["function"]:
                data_grid["function"][ac] = []
            data_grid["function"][ac].append(ann)

    return



def load_ac2pdb(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?pdburi"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac rdfs:seeAlso ?pdburi ."
    qs += " ?pdburi up:database <http://purl.uniprot.org/database/PDB> ."
    qs += " ?ac rdf:type up:Protein ."
    qs += "}"


    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            pdb_id = result["pdburi"]["value"].split("/")[-1]
            if ac not in data_grid["ac2pdb"]:
                data_grid["ac2pdb"][ac] = []
            data_grid["ac2pdb"][ac].append(pdb_id)

    return


def load_ac2reactome(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?reactomeuri ?reactomename"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac rdfs:seeAlso ?reactomeuri ."
    qs += " ?reactomeuri up:database <http://purl.uniprot.org/database/Reactome> ."
    qs += " ?reactomeuri rdfs:comment ?reactomename ."
    qs += " ?ac rdf:type up:Protein ."
    qs += "}"


    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            reactome_id = result["reactomeuri"]["value"].split("/")[-1]
            reactome_label = result["reactomename"]["value"]
            if ac not in data_grid["reactome"]:
                data_grid["reactome"][ac] = []
            data_grid["reactome"][ac].append({"id":reactome_id, "label":reactome_label})

    return


def load_refseq2isoform(sheet_obj, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?resourceuri ?isoformuri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?resourceuri rdfs:seeAlso ?isoformuri ."
    qs += " ?resourceuri rdf:type up:Resource ."
    qs += " ?isoformuri rdf:type up:Simple_Sequence ."
    qs += "}"
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            res_ac = result["resourceuri"]["value"].split("/")[-1]
            res_type = result["resourceuri"]["value"].split("/")[-2]
            if res_type == "refseq":
                isoform = result["isoformuri"]["value"].split("/")[-1]
                if isoform not in sheet_obj:
                    sheet_obj[isoform] = []
                sheet_obj[isoform].append(res_ac)

    return



def load_participants_reactome(sheet_obj, species_list):

    
    seen = {}
    for species in species_list:
        graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
        qs = prefixes
        qs += " SELECT"
        qs += " ?reactionuri ?inputuri ?outputuri ?inputxrefidtype ?outputxrefidtype ?inputxrefid ?outputxrefid ?inputname ?outputname "
        qs += " FROM <%s>" % (graph_uri)
        qs += " WHERE {"
        qs += " ?reactionuri gly:rxnInput ?inputuri ."
        qs += " ?reactionuri gly:rxnOutput ?outputuri ."
        qs += " ?inputuri gly:xrefIdentifier ?inputxrefiduri ."
        qs += " ?outputuri gly:xrefIdentifier ?outputxrefiduri ."
        qs += " ?inputxrefiduri gly:xrefIdType ?inputxrefidtype ."
        qs += " ?outputxrefiduri gly:xrefIdType ?outputxrefidtype ."
        qs += " ?inputxrefiduri gly:xrefId ?inputxrefid ."
        qs += " ?outputxrefiduri gly:xrefId ?outputxrefid ."
        qs += " ?inputuri gly:participantName ?inputname ."
        qs += " ?outputuri gly:participantName ?outputname ."
        qs += " ?reactionuri rdf:type gly:Reaction_Annotation ."
        qs += "}"

        limit = 100000
        total = 0
        for i in xrange(0, 10000):
            newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
            sparql.setQuery(newqs)
            results = sparql.query().convert()
            n = len(results["results"]["bindings"])
            if n == 0:
                break
            total += n
            for result in results["results"]["bindings"]:
                reaction_id = result["reactionuri"]["value"].split("/")[-1]
                input_id = result["inputuri"]["value"].split("/")[-1]
                output_id = result["outputuri"]["value"].split("/")[-1]
                input_xref_id = result["inputxrefid"]["value"].split("/")[-1]
                output_xref_id = result["outputxrefid"]["value"].split("/")[-1]
                input_xref_id_type = result["inputxrefidtype"]["value"]
                output_xref_id_type = result["outputxrefidtype"]["value"]
                input_name = result["inputname"]["value"]
                output_name = result["outputname"]["value"]
                if reaction_id not in sheet_obj["participants"]:
                    sheet_obj["participants"][reaction_id] = {"input":[], "output":[]}
                
                o = {"id":input_id, "xrefid":input_xref_id, "xreftype":input_xref_id_type, "name":input_name}
                o_str = json.dumps(o)
                if o_str not in seen:
                    sheet_obj["participants"][reaction_id]["input"].append(o)
                    seen[o_str] = True
                o = {"id":output_id, "xrefid":output_xref_id, "xreftype":output_xref_id_type, "name":output_name}
                o_str = json.dumps(o)
                if o_str not in seen:
                    sheet_obj["participants"][reaction_id]["output"].append(o)
                    seen[o_str] = True

    return


def load_enzymes_rhea(row_list, species):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?proteinuri ?activityuri ?reactionuri "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?proteinuri up:annotation ?annotationuri ."
    qs += " ?annotationuri up:catalyticActivity ?activityuri ."
    qs += " ?annotationuri rdf:type up:Catalytic_Activity_Annotation ."
    qs += " ?activityuri up:catalyzedReaction ?reactionuri ."
    qs += " ?activityuri rdf:type up:Catalytic_Activity."
    qs += "}"
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            uniprotkb_ac = result["proteinuri"]["value"].split("/")[-1]
            reaction_id = result["reactionuri"]["value"].split("/")[-1]
            row_list.append([uniprotkb_ac, reaction_id])
    return


def load_enzymes_reactome(row_list, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?proteinuri ?reactionuri ?enzymeuri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?proteinuri up:enzyme ?enzymeuri ."
    qs += " ?proteinuri up:annotation ?reactionuri ."
    qs += " ?reactionuri rdf:type gly:Reaction_Annotation ."
    qs += "}"
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            uniprotkb_ac = result["proteinuri"]["value"].split("/")[-1]
            reaction_id = result["reactionuri"]["value"].split("/")[-1]
            row_list.append([uniprotkb_ac, reaction_id])

    return


def load_pathways_reactome(sheet_obj, species):


    pathway2reaction = {}
    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?reactionuri ?reactionname ?pathwayuri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?reactionuri gly:pathway ?pathwayuri ."
    qs += " ?reactionuri rdf:type gly:Reaction_Annotation ."
    qs += "}"
    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            pathway_id = result["pathwayuri"]["value"].split("/")[-1]
            reaction_id = result["reactionuri"]["value"].split("/")[-1]
            if pathway_id not in pathway2reaction:
                pathway2reaction[pathway_id] = []
            if reaction_id not in pathway2reaction[pathway_id]:
                pathway2reaction[pathway_id].append(reaction_id)



    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?pathwayuri ?pathwayname ?pathwaysummary"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?pathwayuri gly:pathwayName ?pathwayname ."
    qs += " ?pathwayuri gly:pathwaySummary ?pathwaysummary ."
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            pathway_id = result["pathwayuri"]["value"].split("/")[-1]
            pathway_name = result["pathwayname"]["value"]
            pathway_summary = result["pathwaysummary"]["value"]
            reactions = ""
            if pathway_id in pathway2reaction:
                reactions = ",".join(pathway2reaction[pathway_id])
            sheet_obj.append([pathway_id, pathway_name, pathway_summary, reactions])

    return


def load_reactions(sheet_obj, species_list):


    reaction2location = {}
    for species in species_list:
        graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
        qs = prefixes
        qs += " SELECT"
        qs += " ?reactionuri ?cellularlocation "
        qs += " FROM <%s>" % (graph_uri)
        qs += " WHERE {"
        qs += " ?reactionuri gly:cellularLocation ?cellularlocation ."
        qs += " ?reactionuri rdf:type gly:Reaction_Annotation ."
        qs += "}"
        limit = 100000
        total = 0
        for i in xrange(0, 10000):
            newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
            sparql.setQuery(newqs)
            results = sparql.query().convert()
            n = len(results["results"]["bindings"])
            if n == 0:
                break
            total += n
            for result in results["results"]["bindings"]:
                reaction_id = result["reactionuri"]["value"].split("/")[-1]
                cellular_location = result["cellularlocation"]["value"]
                if reaction_id not in reaction2location:
                    reaction2location[reaction_id] = []   
                reaction2location[reaction_id].append(cellular_location)

    reaction2pmid = {}
    for species in species_list:
        graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
        qs = prefixes
        qs += " SELECT"
        qs += " ?reactionuri ?citeuri "
        qs += " FROM <%s>" % (graph_uri)
        qs += " WHERE {"
        qs += " ?reactionuri up:citation ?citeuri ."
        qs += " ?reactionuri rdf:type gly:Reaction_Annotation ."
        qs += "}"
        limit = 100000
        total = 0
        for i in xrange(0, 10000):
            newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
            sparql.setQuery(newqs)
            results = sparql.query().convert()
            n = len(results["results"]["bindings"])
            if n == 0:
                break
            total += n
            for result in results["results"]["bindings"]:
                reaction_id = result["reactionuri"]["value"].split("/")[-1]
                pmid = result["citeuri"]["value"].split("/")[-1]
                if reaction_id not in reaction2pmid:
                    reaction2pmid[reaction_id] = []
                reaction2pmid[reaction_id].append(pmid)
    




    seen = {}
    for species in species_list:
        graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
        qs = prefixes
        qs += " SELECT"
        qs += " ?reactionuri ?reactionname ?pathwayuri ?reactionsummary ?reactiondb ?reactionequation"
        qs += " FROM <%s>" % (graph_uri)
        qs += " WHERE {"
        qs += " ?reactionuri rdf:type gly:Reaction_Annotation ."
        qs += " optional { ?reactionuri gly:reactionDatabase ?reactiondb . } "
        qs += " optional { ?reactionuri gly:rxnName ?reactionname . } "
        qs += " optional { ?reactionuri gly:rxnSummary ?reactionsummary . } "
        qs += " optional { ?reactionuri gly:pathway ?pathwayuri . } "
        qs += " optional { ?reactionuri gly:equation ?reactionequation . } "
        qs += "}"

        limit = 100000
        total = 0
        for i in xrange(0, 10000):
            newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
            sparql.setQuery(newqs)
            results = sparql.query().convert()
            n = len(results["results"]["bindings"])
            if n == 0:
                break
            total += n
            for result in results["results"]["bindings"]:
                reaction_id = result["reactionuri"]["value"].split("/")[-1]
                reaction_name, reaction_summary, reaction_db,reaction_equation,pathway_id = "","","","",""
                if "reactiondb" in result:
                    reaction_db = result["reactiondb"]["value"].split("/")[-1]
                if "reactionname" in result:
                    reaction_name = result["reactionname"]["value"]
                if "reactionequation" in result:
                    reaction_equation = result["reactionequation"]["value"]

                if "pathwayuri" in result:
                    pathway_id = result["pathwayuri"]["value"].split("/")[-1]
                if "reactionsummary" in result:
                    reaction_summary = result["reactionsummary"]["value"]
                location = ", ".join(reaction2location[reaction_id]) if reaction_id in reaction2location else ""

                pmid = ", ".join(reaction2pmid[reaction_id]) if reaction_id in reaction2pmid else ""
                combo_id = "%s %s" % (reaction_id, pmid)
                if combo_id not in seen:
                    seen[combo_id] = True
                    o = {"reactionid":reaction_id, "reactionname":reaction_name, 
                            "cellularlocation":location, "pmid":pmid, 
                            "pathwayid":pathway_id, "summary":reaction_summary,
                            "equation":reaction_equation, "source":reaction_db}
                    sheet_obj["reactions"].append(o)

    return


def load_ac2xref(sheet_obj, species, xref_obj):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?uri ?label" if xref_obj["lblflag"] == True else "?ac ?uri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac rdfs:seeAlso ?uri ."
    qs += " ?uri up:database <http://purl.uniprot.org/database/%s> ." % (xref_obj["dbname"])
    qs += " ?uri rdfs:comment ?label ." if xref_obj["lblflag"] == True else ""
    qs += " ?ac rdf:type up:Protein ."
    qs += "}"
    

    if xref_obj["dbname"] == "Reactome":
        qs = prefixes
        qs += " SELECT"
        qs += " ?ac ?uri ?label" if xref_obj["lblflag"] == True else "?ac ?uri"
        qs += " FROM <%s>" % (graph_uri)
        qs += " WHERE {"
        qs += " ?ac up:annotation ?reactionuri ."
        qs += " ?ac rdf:type up:Protein ."
        qs += " ?reactionuri rdf:type gly:Reaction_Annotation ."
        qs += " ?reactionuri gly:pathway ?uri ."
        qs += " ?uri gly:pathwayName ?label ."
        qs += "}"
    elif xref_obj["dbname"] == "Rhea":
        qs = prefixes
        qs += " SELECT"
        qs += " ?ac ?uri ?label "
        qs += " FROM <%s>" % (graph_uri)
        qs += " WHERE {"
        qs += " ?ac up:annotation ?annotationuri ."
        qs += " ?annotationuri up:catalyticActivity ?activityuri ."
        qs += " ?annotationuri rdf:type up:Catalytic_Activity_Annotation ."
        qs += " ?activityuri up:catalyzedReaction ?uri ."
        qs += " ?activityuri rdf:type up:Catalytic_Activity."
        qs += " ?uri gly:equation ?label ."
        qs += "}"


    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            db_id = result["uri"]["value"].split("/")[-1]
            label = result["label"]["value"] if xref_obj["lblflag"] == True else ""
            if ac not in sheet_obj:
                sheet_obj[ac] = []
            sheet_obj[ac].append({"id":db_id, "label":label})

    return







def load_ac2kegg(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?kegguri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac rdfs:seeAlso ?kegguri ."
    qs += " ?kegguri up:database <http://purl.uniprot.org/database/KEGG> ."
    qs += " ?ac rdf:type up:Protein ."
    qs += "}"


    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            kegg_id = result["kegguri"]["value"].split("/")[-1]
            if ac not in data_grid["kegg"]:
                data_grid["kegg"][ac] = []
            data_grid["kegg"][ac].append({"id":kegg_id, "label":""})

def load_protein_existence(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?pe"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:existence ?pe ."
    qs += " ?ac rdf:type up:Protein ."
    qs += "}"

    pe_dict = {
        "Evidence_at_Protein_Level_Existence":1
        ,"Evidence_at_Transcript_Level_Existence":2
        ,"Inferred_from_Homology_Existence":3
        ,"Predicted_Existence":4
        ,"Uncertain_Existence":5
    }

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            pe = result["pe"]["value"].split("/")[-1]
            data_grid["pe"][ac] = pe_dict[pe] if pe in pe_dict else ""
    return



def load_proteinid(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?proteinid "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:mnemonic ?proteinid ."
    qs += " ?ac rdf:type up:Protein ."
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            protein_id = result["proteinid"]["value"].split("/")[-1]
            data_grid["proteinid"][ac] = protein_id

    return


def load_isoformmass(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?isoform ?mass"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?isoform up:mass ?mass ."
    qs += " ?isoform rdf:type up:Simple_Sequence ."
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            isoform = result["isoform"]["value"].split("/")[-1]
            mass = result["mass"]["value"].split("/")[-1]
            data_grid["isoformmass"][isoform] = mass
    
    return


def load_isoformlen(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?isoform ?seq"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?isoform rdf:value ?seq ."
    qs += " ?isoform rdf:type up:Simple_Sequence ."
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            isoform = result["isoform"]["value"].split("/")[-1]
            isoform_seq = result["seq"]["value"].split("/")[-1].strip()
            isoform_len = str(len(isoform_seq))
            data_grid["isoformlen"][isoform] = isoform_len

    return


def load_isoformseq(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?isoform ?seq ?ver"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?isoform rdf:value ?seq ."
    qs += " ?isoform up:version ?ver ."
    qs += " ?isoform rdf:type up:Simple_Sequence ."
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            isoform = result["isoform"]["value"].split("/")[-1]
            isoform_seq = result["seq"]["value"].split("/")[-1].strip()
            isoform_ver = result["ver"]["value"].split("/")[-1].strip()
            data_grid["isoformseq"][isoform] = isoform_seq
            data_grid["isoformver"][isoform] = isoform_ver


    return




def load_test(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?isoformuri ?trsid ?pepid ?chrid ?startpos ?endpos "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?isoformuri gly:ensTranscript ?trsid . "
    qs += "?trsid up:translatedTo ?pepid . "
    qs += "?trsid gly:chromosome ?chrid . "
    qs += "?trsid gly:transcriptRange ?rangeuri . "
    qs += "?rangeuri faldo:begin ?beginuri . "
    qs += "?rangeuri faldo:end ?enduri . "
    qs += "?beginuri faldo:position ?startpos . "
    qs += "?enduri faldo:position ?endpos . "
    qs += "}"

    limit = 100000
    total = 0
    for i in xrange(0, 10000):
        newqs = qs + " LIMIT %s OFFSET %s " % (limit,  i*limit)
        sparql.setQuery(newqs)
        results = sparql.query().convert()
        n = len(results["results"]["bindings"])
        if n == 0:
            break
        total += n
        for result in results["results"]["bindings"]:
            ac = result["ac"]["value"].split("/")[-1]
            gene_uri = result["geneuri"]["value"].split("/")[-1].split("#")[-1]
            gene_name = result["genename"]["value"].split("/")[-1]

    return


