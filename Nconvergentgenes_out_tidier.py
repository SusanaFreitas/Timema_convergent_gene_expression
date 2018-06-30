# Nconvergentgenes_out_tidier.py

import sys
import os
import getopt
import numpy

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_dir_name				  = "NOTHINGSET"
output_filename_base      = "NOTHINGSET"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** Nconvergentgenes_out_tidier.py | written by DJP 25/06/18, Edinburgh ****\n")
		print("\nThis script tidies up output dir Nconvergentgenes_out from 10sp_EdgeR_for_randomised_datasets.R")
		print("\n**** USAGE **** \n")

		print("python Nconvergentgenes_out_tidier.py -i [in_dir_name] -o [output base name]\n\n")
		 

		sys.exit(2)
		
	elif opt in ('-i'):
		in_dir_name = arg
	elif opt in ('-o'):
		output_filename_base = arg
	else:
		print("i dont know")
		sys.exit(2)
		

path = in_dir_name
output_file = open(output_filename_base + "_Nconvergentgenes_out_tidied.csv" , "w")
output_file.write("Perm_N,WB,RT,LG\n")

for path, subdirs, files in os.walk(path):
	for name in files:
		infile = open(os.path.join(path, name))

		perm_N = "N_" + name.split("_N_")[1].rstrip("_")
		
		line_N = 0
		for line in infile:
			line_N = line_N + 1
			if line_N > 1:
				output_file.write(perm_N + "," + line)
				
			
print("\n\nFinished, Holston\n\n")
			
		





