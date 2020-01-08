import os,sys
import json
import csv

__version__="1.0"
__status__ = "Dev"


def load_ann(in_file):

    ann = {}
    with open(in_file, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_frame:
            row_count += 1
            if row_count > 1:
                cls = row[0]
                if cls not in ann:
                    ann[cls] = {
                        "label":row[1], "comment":row[2], "definedby":row[3], "parent":row[4], "seealso":row[5]
                    }
    return ann


###############################
def main():

    
    config_obj = json.loads(open("../../conf/config-1.1.json", "r").read())
    path_obj  =  config_obj[config_obj["server"]]["pathinfo"]
    ns_map =  config_obj["nsmap"]
    uri_map = config_obj["urimap"]


    cls_ann = load_ann("generated/datamodel/class-annotation.csv")
    prd_ann = load_ann("generated/datamodel/predicate-annotation.csv")
    selected_ns_list = ["up", "glycan", "gco", "gly", "ens", "faldo", "xsd", "rdfs", "skos"]
    
    #dcterms
    #IRI
    #owl
    #rdf


    in_file = "generated/datamodel/datamodel.csv"
    seen = {"cls":{}, "prd":{}}
    with open(in_file, "r") as FR:
        for line in FR:
            if line[0] == "#":
                continue
            parts = line.strip().split(",")
            prd = parts[0]
            cls_one, cls_two =  parts[1],  parts[2]
            ns_one, ns_two, ns_three = cls_one.split(":")[0], cls_two.split(":")[0], prd.split(":")[0]
            if ns_one in selected_ns_list:
                seen["cls"][cls_one] = 1
            if ns_two in selected_ns_list:
                seen["cls"][cls_two] = 1
            if ns_three in selected_ns_list:
                if prd not in seen["prd"]:
                    seen["prd"][prd] = []
                seen["prd"][prd].append({"domain":cls_one, "range":cls_two})

    #for prd in seen["prd"]:
    #    if len(seen["prd"][prd]) > 1:
    #        print prd, len(seen["prd"][prd])

    cls_tmpl_dict = {
        'label':'\t\t<rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</rdfs:label>\n'
        ,'subclass':'\t\t<rdfs:subClassOf rdf:resource="%s"/>\n'
        ,'comment':'\t\t<rdfs:comment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</rdfs:comment>\n'
        ,'definedby':'\t\t<rdfs:isDefinedBy rdf:resource="%s"/>\n'
        ,'seealso':'\t\t<rdfs:seeAlso rdf:resource="%s"/>\n'
    }
    
    prd_tmpl_dict = {
        'label':'\t\t<rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</rdfs:label>\n'
        ,'comment':'\t\t<rdfs:comment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</rdfs:comment>\n'
        ,'definedby':'\t\t<rdfs:isDefinedBy rdf:resource="%s"/>\n'
        ,'seealso':'\t\t<rdfs:seeAlso rdf:resource="%s"/>\n'
        ,'domain':'\t\t<rdfs:domain rdf:resource="%s"/>\n'
        ,'range':'\t\t<rdfs:range rdf:resource="%s"/>\n'
    }


    out_file = "generated/datamodel/glygen.owl"
    FW = open(out_file, "w")
    with open("generated/datamodel/header.txt", "r") as FR:
        for line in FR:
            FW.write("%s" % (line))


    for cls in seen["cls"]:
        ns = cls.split(":")[0]
        class_name = cls.split(":")[1]
        cls_ann[cls] = cls_ann[cls] if cls in cls_ann else "N/A"
        url = ns_map[ns] +  class_name
        buff = '\t<owl:Class rdf:about="%s">\n' % (url)
        for k in cls_tmpl_dict:
            buff += cls_tmpl_dict[k] % (cls_ann[cls][k]) if k in cls_ann[cls] else ""
        buff += '\t</owl:Class>\n\n'
        FW.write("%s" % (buff))

    for prd in seen["prd"]:
        for o in seen["prd"][prd]:
            ns = prd.split(":")[0]
            prd_name = prd.split(":")[1]
            url = ns_map[ns] +  prd_name
            if prd in prd_ann:
                for k in prd_tmpl_dict:
                    buff += prd_tmpl_dict[k] % (prd_ann[prd][k]) if k in prd_ann[prd] else ""
            buff = '\t<owl:ObjectProperty rdf:about="%s">\n' % (url)
            domain_url = o["domain"]
            range_url = o["range"]
            if o["domain"].find(":") != -1:
                ns,cls_name = o["domain"].split(":")
                domain_url = ns_map[ns] + cls_name
            if o["range"].find(":") != -1:
                ns,cls_name = o["range"].split(":")
                range_url = ns_map[ns] + cls_name
            buff += prd_tmpl_dict["domain"] % (domain_url)
            buff += prd_tmpl_dict["range"] % (range_url)
            buff += '\t</owl:ObjectProperty>\n\n'
            FW.write("%s" % (buff))




    FW.write("%s" % ("</rdf:RDF>"))
    FW.close()


if __name__ == '__main__':
	main()
