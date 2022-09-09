import os,sys
import json

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


__version__="1.0"
__status__ = "Dev"


###############################
def main():


    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)

    config_json = json.loads(open("../conf/config-1.0.json", "r").read())

    glycandb_file = config_json["files"]["glycandb"]
    glycandb_obj = json.loads(open(glycandb_file, "r").read())

    glycan_uri = "http://purl.glygen.org/glycan"
    motif_uri = "http://purl.glygen.org/motif"
    reaction_uri = "http://purl.glygen.org/reaction"


    #keys
    #[id, iupac, wurcs,glycoct, mass, number_monosaccharides, 
    # classification, species, enzyme, motifs, glycoprotein, crossref]


    for glycan_id in glycandb_obj:
        obj = glycandb_obj[glycan_id]
        #print json.dumps(obj.keys(), indent=4)
        #sys.exit()

        #Type triple
        s = "<%s/%s/>" % (glycan_uri,glycan_id)
        p = "<%s#%s>" % (config_json["namespaces"]["rdf"], "type")
        o = "<%s#%s>" % (config_json["namespaces"]["glycan"], "saccharide")
        nt = "%s %s %s ." % (s, p, o)

        #Sequence triples
        seqformat_list = ["wurcs", "iupac", "glycoct"]
        for seqformat in seqformat_list:
            if seqformat not in obj:
                continue
            if obj[seqformat] == "":
                continue

            p = "<%s#%s>" % (config_json["namespaces"]["glycan"],"has_glycosequence")
            o = "<%s/%s/%s/>" % (glycan_uri,glycan_id,seqformat)
            nt = "%s %s %s ." % (s, p, o)
            print nt
        
            s = o
            p = "<%s#%s>" % (config_json["namespaces"]["rdf"], "type")
            o = "<%s#%s>" % (config_json["namespaces"]["glycan"], "glycosequence")
            nt = "%s %s %s ." % (s, p, o)
            print nt
            
            p = "<%s#%s>" % (config_json["namespaces"]["glycan"], "in_carbohydrate_format")
            o = "<%s#%s>" % (config_json["namespaces"]["glycan"], "carbohydrate_format_"+seqformat)
            nt = "%s %s %s ." % (s, p, o)
            print nt

            p = "<%s#%s>" % (config_json["namespaces"]["glycan"], "has_sequence")
            o = "\"%s\"" % (obj[seqformat])
            nt = "%s %s %s ." % (s, p, o)
            print nt
    

        #Enzyme triples
        for eobj in obj["enzyme"]:
            s = "<%s/%s/>" % (glycan_uri,glycan_id)
            p = "<%s#%s>" % (config_json["namespaces"]["glycan"], "synthesized_by")
            o = "<%s/%s/>" % (reaction_uri,"xxx_instanceof_reaction")
            nt = "%s %s %s ." % (s, p, o)
            print "flag-1", nt


        #Motifs triples
        for mobj in obj["motifs"]:
            s = "<%s/%s/>" % (glycan_uri,glycan_id)
            p = "<%s#%s>" % (config_json["namespaces"]["glycan"], "has_motif")
            o = "<%s/%s/>" % (motif_uri,mobj["id"])
            nt = "%s %s %s ." % (s, p, o)
           
            #print "flag-1", nt
            s = o
            p = "<%s#%s>" % (config_json["namespaces"]["rdf"], "type")
            o = "<%s#%s>" % (config_json["namespaces"]["glycan"], "glycan_motif")
            nt = "%s %s %s ." % (s, p, o)
            #print "flag-2", nt

            p = "<%s#%s>" % (config_json["namespaces"]["glycan"], "has_glycoconjugate_sequence")
            o = "<%s#%s>" % (config_json["namespaces"]["glycan"], "xxxx_instanceof_glycoconjugate_sequence")
            nt = "%s %s %s ." % (s, p, o)
            #print "flag-3", nt
                        
            p = "<%s#%s>" % (config_json["namespaces"]["glycan"], "has_glycosequence")
            o = "<%s#%s>" % (config_json["namespaces"]["glycan"], "xxxx_instanceof_glycosequence")            
            nt = "%s %s %s ." % (s, p, o)
            #print "flag-4", nt


            #glycan:has_glycoconjugate_sequence  glycan:glycoconjugate_sequence
            #glycan:has_glycosequence    glycan:glycosequence
            #glycan:contained_in glycan:saccharide

                            


    




if __name__ == '__main__':
	main()
