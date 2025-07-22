import os,sys
import string
import csv
import json
from optparse import OptionParser
import numpy


__version__="1.0"
__status__ = "Dev"




###############################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-c", "--clusterid", action = "store", dest = "clusterid", help = "Cluster ID <1, 2, 3>")


    (options,args) = parser.parse_args()
    for file in ([options.clusterid]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    cluster_id = options.clusterid

    in_file_one = "tmp/entrotypes.csv"
    #in_file_two = "tmp/MetaHIT_SangerSamples.genus.txt"
    in_file_two = "tmp/MetaHIT_SangerSamples.genus.filtered.csv"


    group_dict = {}
    with open(in_file_one, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        rowcount = 0
        for row in data_frame:
            rowcount += 1
            if rowcount > 1:
                if row[1] not in group_dict:
                    group_dict[row[1]] = []
                group_dict[row[1]].append(row[0].strip())

    sample_list = []
    sheet_obj = {}
    with open(in_file_two, 'r') as FR:
        data_frame = csv.reader(FR, delimiter=',', quotechar='"')
        rowcount = 0
        for row in data_frame:
            rowcount += 1
            if rowcount == 1:
                sample_list = row
            elif rowcount > 2:
                for j in xrange(1,len(row)):
                    genus = row[0]
                    sample = sample_list[j].strip().replace("-", ".")
                    if genus not in sheet_obj:
                        sheet_obj[genus] = {}
                    ab = float(row[j])
                    sheet_obj[genus][sample] = ab


    sample_sum = {}
    for genus in sheet_obj:
        for sample in sheet_obj[genus]:
            if sample not in sample_sum:
                sample_sum[sample] = 0.0
            sample_sum[sample] += sheet_obj[genus][sample]



    grp_dict = {}
    for genus in sheet_obj:
        for grp in group_dict:
            ab_list = []
            for sample in group_dict[grp]:
                ab = float(sheet_obj[genus][sample]) if sample in sheet_obj[genus] else 0.0
                ab_list.append(ab)
                #ab_list.append(ab/sample_sum[sample])

            if grp not in grp_dict:
                grp_dict[grp] = {}
            grp_dict[grp][genus] = numpy.var(ab_list)


    from collections import OrderedDict
    sorted_dict = OrderedDict(sorted(grp_dict[cluster_id].items(), key=lambda x: x[1]))

    for genus in sorted_dict:
        print sorted_dict[genus], genus



if __name__ == '__main__':
    main()








