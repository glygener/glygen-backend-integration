import os,sys
import util
import json

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


__version__="1.0"
__status__ = "Dev"


###############################
def main():



	usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")
        parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

        (options,args) = parser.parse_args()
        for file in ([options.configfile, options.species]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

        configJson = json.loads(open(options.configfile, "r").read())
        speciesName = options.species

        queryJson = configJson["queries_2"]


	orgName = configJson["orgname"][speciesName]
	inputFile = configJson["files"][speciesName]["idmapfile"]
	outFile1 = configJson["files"][speciesName]["fasta_canonical"]
	logFile1 = configJson["files"][speciesName]["log_canonical"]
	outFile2 = configJson["files"][speciesName]["fasta_all"]
        logFile2 = configJson["files"][speciesName]["log_all"]
	queryJson = configJson["queries_2"]
	

	seen = {}
	seen["desc"] = {}
        queryString = queryJson["desc"]
        sparql = SPARQLWrapper("http://localhost:8890/sparql")
        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        for entry in res['results']['bindings']:
                ac = entry["s"]["value"].split("/")[-1].split("#")[0]
                desc = entry["o"]["value"]
                seen["desc"][ac] = desc

	seen["sv"] = {}
        queryString = queryJson["sv"]
        sparql = SPARQLWrapper("http://localhost:8890/sparql")
        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        for entry in res['results']['bindings']:
                ac = entry["s"]["value"].split("/")[-1]
                sv = entry["o"]["value"]
                seen["sv"][ac] = sv


	seen["pe"] = {}
        queryString = queryJson["pe"]
        sparql = SPARQLWrapper("http://localhost:8890/sparql")
        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        for entry in res['results']['bindings']:
                ac = entry["s"]["value"].split("/")[-1]
                pe = entry["o"]["value"].split("/")[-1]
                seen["pe"][ac] = 1 if pe == "Evidence_at_Protein_Level_Existence" else 0

	
	seen["seqid"] = {}
        queryString = queryJson["seqid"]
        sparql = SPARQLWrapper("http://localhost:8890/sparql")
        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        for entry in res['results']['bindings']:
                ac = entry["s"]["value"].split("/")[-1]
                seqid = entry["o"]["value"]
                seen["seqid"][ac] = seqid

	seen["seq"] = {}
	seen["has_sequence"] = {}
        queryString = queryJson["seq"]
        sparql = SPARQLWrapper("http://localhost:8890/sparql")
        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        for entry in res['results']['bindings']:
                isoform = entry["s"]["value"].split("/")[-1].split("#")[0]
                seq = entry["o"]["value"]
                seen["seq"][isoform] = seq
			

        lcount = 0
	seen["canonical"] = {} 
        with open(inputFile, "r") as FR:
                for line in FR:
                        lcount += 1
                        if lcount == 1:
                                continue
                        canonical, status, gene, isoforms1, isoforms2 = line.strip().split(",")
			if canonical not in seen["canonical"]:
				seen["canonical"][canonical] = {}
			seen["canonical"][canonical]["status"] = status
			seen["canonical"][canonical]["sprotisoforms"] = isoforms1  
			seen["canonical"][canonical]["tremblisoforms"] = isoforms2  
			seen["canonical"][canonical]["gene"] = gene

	FW1 = open(outFile1, "w")
	FW2 = open(logFile1, "w")
	seen["done"] = {}
	for canonical  in seen["canonical"]:
		if canonical in seen["done"]:
			continue
		ac = canonical.split("-")[0]
		seqid = seen["seqid"][ac] if ac in seen["seqid"] else "-"
		desc = seen["desc"][ac] if ac in seen["desc"] else "-"
		pe = seen["pe"][ac] if ac in seen["pe"] else "-"
		sv = seen["sv"][ac] if ac in seen["sv"] else "-"
		if canonical not in  seen["seq"]:
			FW2.write("no-sequence exists for canonical=%s\n" %(canonical))
			continue
		seq = seen["seq"][canonical]
		acLbl = ""
		if seen["canonical"][canonical]["status"] == "reviewed":
			acLbl = "sp|%s|%s" % (canonical,seqid)
		else:
			acLbl = "tr|%s|%s" % (canonical,seqid)
		descLbl = "%s OS=%s GN=%s PE=%s SV=%s" %(desc, orgName, gene, pe, sv)
		seqObj = SeqRecord(Seq(seq,IUPAC.protein),id=acLbl, 
					name=acLbl, description=descLbl)
		FW1.write("%s\n" % (seqObj.format("fasta")))
		seen["done"][canonical] = True
	FW1.close()
	FW2.close()

	FW1 = open(outFile2, "w")
        FW2 = open(logFile2, "w")
	seen["done"] = {}
        for canonical  in seen["canonical"]:
                for isoform in seen["canonical"][canonical]["sprotisoforms"].split("|"):
			if isoform in seen["done"] or isoform == "":
                        	continue
			ac = isoform.split("-")[0]
                	seqid = seen["seqid"][ac] if ac in seen["seqid"] else "-"
                	desc = seen["desc"][ac] if ac in seen["desc"] else "-"
                	pe = seen["pe"][ac] if ac in seen["pe"] else "-"
                	sv = seen["sv"][ac] if ac in seen["sv"] else "-"
                	if isoform not in  seen["seq"]:
				FW2.write("no-sequence exists for isoform=%s\n" % (isoform))
                		continue
			seq = seen["seq"][isoform]
			acLbl = "sp|%s|%s|%s" % (isoform,canonical,seqid)
                	descLbl = "%s OS=%s GN=%s PE=%s SV=%s" %(desc, orgName, gene, pe, sv)
                	seqObj = SeqRecord(Seq(seq,IUPAC.protein),id=acLbl,
                        	                name=acLbl, description=descLbl)
                	FW1.write("%s\n" % (seqObj.format("fasta")))
                	seen["done"][isoform] = True

		for isoform in seen["canonical"][canonical]["tremblisoforms"].split("|"):
                        if isoform in seen["done"] or isoform == "":
                                continue
                        ac = isoform.split("-")[0]
                        seqid = seen["seqid"][ac] if ac in seen["seqid"] else "-"
                        desc = seen["desc"][ac] if ac in seen["desc"] else "-"
                        pe = seen["pe"][ac] if ac in seen["pe"] else "-"
                        sv = seen["sv"][ac] if ac in seen["sv"] else "-"
                        if isoform not in  seen["seq"]:
                                FW2.write("no-sequence exists for TREMBL isoform=%s\n" %(isoform))
                        	continue	
			seq = seen["seq"][isoform]
                        acLbl = "tr|%s|%s|%s" % (isoform,canonical,seqid)
                        descLbl = "%s OS=%s GN=%s PE=%s SV=%s" %(desc, orgName, gene, pe, sv)
                        seqObj = SeqRecord(Seq(seq,IUPAC.protein),id=acLbl,
                                                name=acLbl, description=descLbl)
                        FW1.write("%s\n" % (seqObj.format("fasta")))
                        seen["done"][isoform] = True
	FW1.close()
	FW2.close()

		#Cannonical
		#>sp|P04278|SHBG_HUMAN Sex hormone-binding globulin OS=Homo sapiens GN=SHBG PE=1 SV=2 
		#SP isoform
		#>sp|P04278-2|SHBG_HUMAN Isoform 2 of Sex hormone-binding globulin OS=Homo sapiens GN=SHBG
		#Tr isoform (where M12345 is a Tr isoform mapped to cannonical SP)
		#>tr|P04278-M12345|SHBG_HUMAN isoform GN=SHBG
		#Tr isoform (where G56789 is a Tr isoform mapped to cannonical Tr T12345)
		#>tr|T12345-G56789|T12345_HUMAN isoform GN=SHBG
 






if __name__ == '__main__':
	main()


