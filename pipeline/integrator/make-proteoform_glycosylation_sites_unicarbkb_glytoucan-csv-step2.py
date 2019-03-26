import os,sys
import json


from optparse import OptionParser
from SPARQLWrapper import SPARQLWrapper, JSON 


__version__="1.1"
__status__ = "Dev"


"""
BCF Nov. 20, 2018
Takes list file of motifs associated with unicarbkb ids and creates a list of N-Glycans and O-Glycans for each unicarbkb id
This script supersedes steps 2a and 2b since we now recieve this informtion directly from unicarbkb and do not need to do
web scraping to find it.
"""
##############################
def addType(dictionary, uckb_id, motif):
    #add to dictionary if not there
    if not(uckb_id in dictionary):
        dictionary[uckb_id] = {"N-Glycan": False, "O-Glycan": False}
    if "N-Glycan" in motif:
        dictionary[uckb_id]["N-Glycan"] = True
    if "O-Glycan" in motif:
        dictionary[uckb_id]["O-Glycan"] = True
        

###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--configfile",action="store",dest="configfile",help="config JSON file")
    parser.add_option("-s","--species",action="store",dest="species",help="human/mouse")

    (options,args) = parser.parse_args()
    for file in ([options.species, options.configfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    config_obj = json.loads(open(options.configfile, "r").read())
    species = options.species
    tax_id = "9606" if species == "human" else "10090" #TODO add support for other species


    in_file = config_obj["pathinfo"]["downloads"] + "/unicarbkb/%s_motif_current.txt" % (species)
#in_file format
### README
#Human glycans with GlyToucan Ids registed in UniCarbKB (14th Sept 2018)
#Class processed internally and validated with GlyToucan motifs
#Note duplication of ids. This is being checked (see duplicated toucanIds.txt)
#
#10 , G80903UK , O-Glycan core 6 Fuzzy
#10 , G80903UK , O-Glycan core 1 fuzzy
#10 , G80903UK , O-Glycan core 2
    uckb_ids = {} #dict will store tuple of bools for glycosylation_types (N-Glycan or O-Glycan)
    with open(in_file, 'r') as inpt:
        for line in inpt:
            spl = line.split(",")
            if ((len(spl) == 3) and (spl[0].rstrip().isdigit())): #3 columns for now, uckb_id is integer but file has extra white space at end
                if ("-Glycan" in spl[2]):
                    addType(uckb_ids, spl[0], spl[2]) 
            
            


    out_file1 = config_obj["pathinfo"]["intermediate"] + "/%s_glycosylation_types.csv" % (species)
    out_file2 = config_obj["pathinfo"]["intermediate"] + "/%s_glycosylation_types.log" % (species)
    FWcsv = open(out_file1, "w")
    FWlog = open(out_file2, "w")

    count_no_type = 0
    count_both_type = 0
   
    row = ["uckb_id", "glycosylation_type"]
    FWcsv.write("\"%s\"\n" % ("\",\"".join(row)))
    #FWlog.write("\"%s\"\n" % ("\",\"".join(row)))
    for uckb_id, glyTypes in uckb_ids.iteritems():
        types_this_glycan = 0
        for glycan_type in ["N-Glycan", "O-Glycan"]: 
            if glyTypes[glycan_type]:
                types_this_glycan += 1
                row = [uckb_id, glycan_type]
                FWcsv.write("\"%s\"\n" % ("\",\"".join(row)))
                #FWlog.write("\"%s\"\n" % ("\",\"".join(row)))
        if (types_this_glycan == 2):   #both types, write to csv and log, filtered later
            count_both_type += 1
            row = [uckb_id, "both types"]
            FWlog.write("\"%s\"\n" % ("\",\"".join(row)))        
        elif (types_this_glycan == 0):   #no type, write to log only
            count_no_type += 1
            row = [uckb_id, "no types"]
            FWlog.write("\"%s\"\n" % ("\",\"".join(row)))        

    FWcsv.close()
    FWlog.close()

    print "%d glycans with no type \n %d glycans with both types \n see %s for more info" % (count_no_type, count_both_type, out_file2)




if __name__ == '__main__':
	main()
