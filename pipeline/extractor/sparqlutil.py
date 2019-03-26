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
    PREFIX gly: <http://glygen-vm-prd.biochemistry.gwu.edu/ontology/>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX oboinowl: <http://www.geneontology.org/formats/oboInOwl#>
"""

def load_do_mapping(data_grid, species):


    graph_uri = "http://sparql.glygen.org#disease"
    qs = prefixes
    qs += " SELECT"
    qs += "?douri ?doval ?xrefval ?doname ?doaltname ?dodef "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?douri oboinowl:id ?doval . "
    qs += "?douri oboinowl:hasDbXref ?xrefval . "
    qs += "?douri rdfs:label ?doname . "
    qs += "?douri oboinowl:hasExactSynonym ?doaltname . "
    qs += "?douri obo:IAO_0000115 ?dodef . "
    qs += "}"

    xobj_list_one = [
        {"prefix":"ICD10CM:", "mapname":"doid2icd10cmid"},
        {"prefix":"ICD9CM:", "mapname":"doid2icd9cmid"},
        {"prefix":"KEGG:", "mapname":"doid2keggid"},
        {"prefix":"MESH:", "mapname":"doid2meshid"},
        {"prefix":"UMLS_CUI:", "mapname":"doid2umlsid"}
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
            xref_val = result["xrefval"]["value"].strip()
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

            do_name = result["doname"]["value"].strip()
            if do_id not in data_grid["doid2name"]:
                data_grid["doid2name"][do_id] = []
            data_grid["doid2name"][do_id].append(do_name)

            do_altname = result["doaltname"]["value"].strip()
            if do_id not in data_grid["doid2altname"]:
                data_grid["doid2altname"][do_id] = []
            data_grid["doid2altname"][do_id].append(do_altname)

            do_def = result["dodef"]["value"].strip()
            if do_id not in data_grid["doid2def"]:
                data_grid["doid2def"][do_id] = []
            data_grid["doid2def"][do_id].append(do_def)


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

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?citeuri ?journalname ?journaltitle "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += "?ac up:citation ?citeuri . "
    qs += "?citeuri up:name ?journalname . "
    qs += "?citeuri up:title ?journaltitle . "
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
                journal_title = result["journaltitle"]["value"].strip()
                if ac not in data_grid["citelist"]:
                    data_grid["citelist"][ac] = []
                o = {"pmid":pm_id, "journalname":journal_name, "journaltitle":journal_title}
                data_grid["citelist"][ac].append(o)
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


def load_isoforminfo(data_grid, species):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?isoform ?status ?iscanon"
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
            isoform = result["isoformuri"]["value"].split("/")[-1]
            if isoform not in data_grid["locusinfo"]:
                data_grid["locusinfo"][isoform] = {
                    "chrid":result["chrid"]["value"].strip()
                    ,"trsid":result["trsid"]["value"].split("/")[-1] 
                    ,"pepid":result["pepid"]["value"].split("/")[-1]
                    ,"strand":"x"
                    ,"startpos":result["startpos"]["value"]
                    ,"endpos":result["endpos"]["value"]
                }

    return



def load_recnames(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?fullname"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:recommendedName ?nameuri ."
    qs += " ?nameuri up:fullName ?fullname ."
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
            full_name = result["fullname"]["value"]
            data_grid["recnames"][ac] = {"fullname":full_name}
  

    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?shortname"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:recommendedName ?nameuri ."
    qs += " ?nameuri up:shortName ?shortname ."
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
            short_name = result["shortname"]["value"]
            if ac in data_grid["recnames"]:
                data_grid["recnames"][ac]["shortname"] = short_name
            else:
                data_grid["recnames"][ac] = {"fullname":"", "shortname":short_name}

    return



def load_altnames(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?fullname"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:alternativeName ?nameuri ."
    qs += " ?nameuri up:fullName ?fullname ."
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
            full_name = result["fullname"]["value"]
            data_grid["altnames"][ac] = {"fullname":full_name}


    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?shortname"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:alternativeName ?nameuri ."
    qs += " ?nameuri up:shortName ?shortname ."
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
            short_name = result["shortname"]["value"]
            if ac in data_grid["altnames"]:
                data_grid["altnames"][ac]["shortname"] = short_name
            else:
                data_grid["altnames"][ac] = {"fullname":"", "shortname":short_name}

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
    qs += " ?annuri rdfs:comment ?comment . "
    qs += " ?annuri rdf:type up:Glycosylation_Annotation . "
    qs += " ?annuri up:attribution ?atturi . "
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
            car_id = result["annuri"]["value"].split("/")[-1]
            car_id = car_id if car_id[0:4] == "CAR_" else ""
            ac = result["ac"]["value"].split("/")[-1]
            comment = result["comment"]["value"]
            evidence = result["evidenceuri"]["value"].split("/")[-1]
            gly_type = comment.strip().split(" ")[0]
            saccharide = comment.strip().split(" ")[1]
            saccharide = saccharide.replace("(", "").replace(")", "")
            start_pos = result["startpos"]["value"]
            end_pos = result["endpos"]["value"]
            source_one = result["citationuri"]["value"].split("/")[-2]
            source_one = source_dict[source_one] if source_one in source_dict else source_one
            source_two = result["citationuri"]["value"].split("/")[-1]
            for aa in aa_dict:
                if comment.find(aa.lower()) != -1:
                    aa = aa_dict[aa]
                    row = [ac, gly_type, saccharide, aa, start_pos]
                    combo_id = ",".join(row)
                    if combo_id not in ev_dict:
                        ev_dict[combo_id] = []
                    ev_dict[combo_id].append([source_one, source_two,evidence, car_id])




    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += " ?ac ?comment ?annuri ?startpos ?endpos "
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac up:annotation ?annuri . "
    qs += " ?annuri rdfs:comment ?comment . "
    qs += " ?annuri rdf:type up:Glycosylation_Annotation . "
    qs += " ?annuri up:range ?rangeuri . "
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
            car_id = result["annuri"]["value"].split("/")[-1]
            car_id = car_id if car_id[0:4] == "CAR_" else ""
            comment = result["comment"]["value"]
            gly_type = comment.strip().split(" ")[0]
            saccharide = comment.strip().split(" ")[1]
            saccharide = saccharide.replace("(", "").replace(")", "")
            start_pos = result["startpos"]["value"]
            end_pos = result["endpos"]["value"]
            for aa in aa_dict:
                if comment.find(aa.lower()) != -1:
                    aa = aa_dict[aa]
                    row_one = [ac, gly_type, saccharide, aa, start_pos]
                    combo_id = ",".join(row_one)
                    if combo_id in ev_dict:
                        for row_two in ev_dict[combo_id]:
                            data_grid["glycosylation"].append(row_one + row_two)
                    else:
                        data_grid["glycosylation"].append(row_one + ["", "", car_id, ""])

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
            data_grid["function"][ac] = ann

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




def load_ac2xref(sheet_obj, species, xref_obj):


    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?uri ?label" if xref_obj["lblflag"] == True else "?ac ?uri"
    qs += " FROM <%s>" % (graph_uri)
    qs += " WHERE {"
    qs += " ?ac rdfs:seeAlso ?uri ."
    qs += " ?uri up:database <http://purl.uniprot.org/database/%s> ." % (xref_obj["dbname"])
    qs += " ?uri rdfs:comment ?label ." if xref_obj["lblflag"] == True else ""
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



def load_proteinid(data_grid, species):

    graph_uri = "http://sparql.glygen.org#uniprot_%s" % (species)
    qs = prefixes
    qs += " SELECT"
    qs += "?ac ?proteinid"
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
            data_grid["isoformseq"][isoform] = isoform_seq


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
            print gene_name, ac, gene_uri

    return


