nsmap = {
    "up":"http://purl.uniprot.org/core/"
    ,"skos":"http://www.w3.org/2004/02/skos/core#"
    ,"rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    ,"owl":"http://www.w3.org/2002/07/owl#"
    ,"ms":"http://www.berkeleybop.org/ontologies/ms.owl"
    ,"xsd":"http://www.w3.org/2001/XMLSchema#"
    ,"rdfs":"http://www.w3.org/2000/01/rdf-schema#"
    ,"glycan":"http://purl.jp/bio/12/glyco/glycan#"
    ,"uk":"http://unicarbkb.org/unicarbkbprotein/"
    ,"faldo":"http://biohackathon.org/resource/faldo#"
    ,"pato":"http://purl.obolibrary.org/obo/uo.owl"
    ,"dcterms":"http://purl.org/dc/terms/"
    ,"bibo":"http://purl.org/ontology/bibo/"
    ,"unicorn":"http://purl.jp/bio/12/glyco/unicorn/"
    ,"foaf":"http://xmlns.com/foaf/0.1/"
    ,"dc":"http://purl.org/dc/elements/1.1/"
    ,"gly":"http://glygen-vm-prd.biochemistry.gwu.edu/ontology/"
}

seen = {"modelprd":{}, "modelcls":{}, "dataprd":{}, "datacls":{}}
model_file = "/var/www/html/datasetviewer/content/datamodel.2.csv"
with open(model_file, "r") as FR:
    for line in FR:
        prd = line.split(",")[0]
        seen["modelprd"][prd] = True
        cls = line.split(",")[1]
        seen["modelcls"][cls] = True
        cls = line.split(",")[2]
        seen["modelcls"][cls] = True



data_file = "downloads/ebi/2018_11/uniprot-proteome-9606-homo-sapiens-2018-11.nt"
with open(data_file, "r") as FR:
    for line in FR:
        if line[0] == "_":
            continue
        
        prd = line.split(" ")[1][1:-1]
        if prd not in seen["dataprd"]:
            newprd = prd
            for prefix in nsmap:
                if prd.find(nsmap[prefix]) != -1:
                    newprd = prefix + ":" + prd.split("/")[-1].split("#")[-1]
                    break
            seen["dataprd"][prd] = True
            if newprd not in seen["modelprd"]:
                print "prd", newprd
   
        if prd  == "http://www.w3.org/1999/02/22-rdf-syntax-ns#type":
            cls = line.split(" ")[2][1:-1]
            if cls not in seen["datacls"]:
                newcls = cls
                for prefix in nsmap:
                    if cls.find(nsmap[prefix]) != -1:
                        newcls = prefix + ":" + cls.split("/")[-1].split("#")[-1]
                        break
                seen["datacls"][cls] = True
                if newcls not in seen["modelcls"]:
                    print "cls", newcls

