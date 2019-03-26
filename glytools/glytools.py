import os,sys
import json
import libgly

from optparse import OptionParser





def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog ")
    parser.add_option("-w","--workflowfile",action="store",dest="workflowfile",help="Workflow JSON file")

    (options,args) = parser.parse_args()
    for file in ([options.workflowfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    wf_json = json.loads(open(options.workflowfile, "r").read())

    workbook = {"sheets":{}}
    
    
    #workbook["sheets"]["sheet_one"] = {}
    #in_file = "workflows/toy.csv"
    #libgly.load_sheet(workbook["sheets"]["sheet_one"], in_file, ",")
    

    libgly.load_workbook(workbook, wf_json["inputfiles"], ",")


    for obj in wf_json["inputfiles"]:
        for file_name in obj["filenamelist"]:
            #print workbook["sheets"][file_name]["data"][0:3]
            print file_name, workbook["sheets"][file_name]["fields"]


if __name__ == '__main__':
    main()



