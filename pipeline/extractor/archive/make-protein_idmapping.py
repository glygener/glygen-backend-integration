import os,sys
import util
import json

from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


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


	inputFile = configJson["files"][speciesName]["groupingfile"]
	outputFile = configJson["files"][speciesName]["idmapfile"]
	queryJson = configJson["queries_1"] 

        tax_id = configJson["taxid"][speciesName]
        sparql = SPARQLWrapper("http://localhost:8890/sparql")
                
	seen = {"sprotac":{}, "tremblac":{}, "gene":{}, "genefromgraph":{}, "isoformlist":{}}
	with open(inputFile, "r") as FR:
		for line in FR:
			parts = line.strip().split("\t")
			ac, gene = parts[0], parts[-1]
			seen["gene"][ac] = gene



	queryString = queryJson["reviewed"]
	queryString = queryString.replace("TAXID", str(tax_id))
	sparql.setQuery(queryString)
	sparql.setReturnFormat(JSON)
	res = sparql.query().convert()
	for entry in res['results']['bindings']:
		ac = entry["s"]["value"].split("/")[-1]
		seen["sprotac"][ac] = 1
        
	queryString = queryJson["unreviewed"]
        queryString = queryString.replace("TAXID", str(tax_id))
        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        for entry in res['results']['bindings']:
                ac = entry["s"]["value"].split("/")[-1]
                seen["tremblac"][ac] = 1
        
	seen["ac2canonical"] = {}
	queryString = queryJson["isoformlist"]
        queryString = queryString.replace("TAXID", str(tax_id))
        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
	for entry in res['results']['bindings']:
                ac = entry["s"]["value"].split("/")[-1]
		isoform = entry["o"]["value"].split("/")[-1]
                if "#" not in ac:
			if ac not in seen["isoformlist"]:
				seen["isoformlist"][ac] = []
			if isoform.find(ac) != -1 and ac not in seen["ac2canonical"]:
				seen["ac2canonical"][ac] = isoform  #defaut canonical to first isoform
                        seen["isoformlist"][ac].append(isoform)



	queryString = queryJson["ac2canonical"]
        queryString = queryString.replace("TAXID", str(tax_id))

        sparql.setQuery(queryString)
        sparql.setReturnFormat(JSON)
        res = sparql.query().convert()
        for entry in res['results']['bindings']:
                isoform = entry["s"]["value"].split("/")[-1]
                canonical = entry["o"]["value"].split("/")[-1]
		ac = isoform.split("-")[0]
		seen["ac2canonical"][ac] = canonical #over write default canonical

	n1, n2 = 0, 0
	seen["gene2canonical"] = {}
	seen["reviewedcanonical"] = {}
	output = {}
	for ac in seen["sprotac"]:
		gene = seen["gene"][ac]
		canonical = seen["ac2canonical"][ac]
		if canonical not in output:
			output[canonical] = {}
		if "sprotsoformlist" not in output[canonical]:
                        output[canonical]["sprotisoformlist"] = []
		output[canonical]["sprotisoformlist"] += seen["isoformlist"][ac]
		output[canonical]["status"] = "reviewed"
		output[canonical]["gene"] = gene
		seen["gene2canonical"][gene] = canonical
		seen["reviewedcanonical"][canonical] = True

	n2a, n2b = 0, 0
	for ac in seen["tremblac"]:
		gene = seen["gene"][ac]
		if gene in seen["gene2canonical"]:
			canonical = seen["gene2canonical"][gene]
		else:
			canonical = seen["isoformlist"][ac][0]
			if canonical not in output:
				output[canonical] = {}
			output[canonical]["status"] = "unreviewed"
		if "tremblisoformlist" not in output[canonical]:
			output[canonical]["tremblisoformlist"] = []
		output[canonical]["tremblisoformlist"] += seen["isoformlist"][ac]
		output[canonical]["gene"] = gene


 
	FW = open(outputFile, "w")
	FW.write("%s,%s,%s,%s,%s\n" % ("canonical_accession", "status", "gene_name","reviewed_isoforms", "unreviewed_isoforms"))
	for canonical in output:
		isoformList1, isoformList2 = [], []
		isoforms1, isoforms2 = "", ""
		if "sprotisoformlist" in output[canonical]:
			isoformList1 = sorted(set(output[canonical]["sprotisoformlist"]))
			isoforms1 = "|".join(isoformList1)
		if "tremblisoformlist" in output[canonical]:
			isoformList2 = sorted(set(output[canonical]["tremblisoformlist"]))
			isoforms2 = "|".join(isoformList2)
		FW.write("%s,%s,%s,%s,%s\n" % (canonical, output[canonical]["status"],
				output[canonical]["gene"], isoforms1, isoforms2))
	FW.close()






if __name__ == '__main__':
	main()


