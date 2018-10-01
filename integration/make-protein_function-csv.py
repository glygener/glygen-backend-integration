import os,sys
import json
import csv
from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON

__version__="1.0"
__status__ = "Dev"

###############################
def main():

    config_json = json.loads(open("../conf/config-1.json", "r").read())
    out_file = config_json["pathinfo"]["unreviewed"]+"/"+"acc_recom_all.csv"
    queries = {
        "q_1":"SELECT ?s, ?o WHERE { ?s ?x <http://purl.uniprot.org/core/Function_Annotation> . ?s <http://www.w3.org/2000/01/rdf-schema#comment> ?o . }",
        "q_2": "SELECT ?s, ?o WHERE { ?s <http://purl.uniprot.org/core/recommendedName> ?x . ?x <http://purl.uniprot.org/core/fullName> ?o . }"
    }

    sparql = SPARQLWrapper("http://localhost:8890/sparql")
    sparql.setQuery(queries["q_2"])
    sparql.setReturnFormat(JSON)
    res = sparql.query().convert()
    acc_name = {}
    uniprotId_list = []
    for entry in res['results']['bindings']:
        subject_id = entry["s"]["value"].split("/")[-1]
        object_id = entry["o"]["value"].split("/")[-1]
        if '#' in subject_id:
            continue
        else:
            acc_name[subject_id] = object_id
    uniprotId_list = uniprotId_list + acc_name.keys()

    sparql.setQuery(queries["q_1"])
    sparql.setReturnFormat(JSON)
    res = sparql.query().convert()
    acc_function = {}
    for entry in res['results']['bindings']:
        subject_id = entry["s"]["value"].split("/")[-1].split("#")[0]
        object_id = entry["o"]["value"]
        if '#' in subject_id:
            continue
        else:
            acc_function[subject_id] = object_id.replace('"', "'")
    uniprotId_list = uniprotId_list + acc_function.keys()
    uniprotId_list = list(set(uniprotId_list))

    headerList = ["uniprotkb_acc_canonical","protein_name","protein_function"]
    with open(outFile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headerList)
        for Id in uniprotId_list:
            row = [Id]
            if Id in acc_name:
                row.append('"'+acc_name[Id]+'"')
            else:
                row.append('')
            if Id in acc_function:
                row.append('"'+acc_function[Id]+'"')
            else:
                row.append('')
            writer.writerow(row)

if __name__ == '__main__':
    main()
