import os,sys
from optparse import OptionParser
import commands



__version__="1.0"
__status__ = "Dev"


###############################
def main():


	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i","--infile",action="store",dest="infile",help="Input file (.nt)")
        parser.add_option("-o","--outdir",action="store",dest="outdir",help="Output folder")

        (options,args) = parser.parse_args()

	for file in ([options.infile, options.outdir]):
		if not (file):
			parser.print_help()
			sys.exit(0)
        nt_file = options.infile
        out_dir = options.outdir
        

        nt_file_name = os.path.basename(nt_file)
	prefix = out_dir + "/" + ('.').join(nt_file_name.split('.')[:-1])
        
        cmd = "rm %s.*.nt" % (prefix)
        x = commands.getoutput(cmd)

        with open(nt_file, "r") as FR:
		i = 1
		j = 1
		out_file = prefix + "." + str(j) + ".nt"
		FW = open(out_file, "w")
		for line in FR:
			FW.write("%s" % (line))
                	if i%5000000 == 0:
				j += 1
                        	FW.close()
				out_file = prefix + "." + str(j) + ".nt"
                		FW = open(out_file, "w")
				print "parsed %s statements" % (i)
			i += 1
		FW.close()


if __name__ == '__main__':
	main()


