import sys
import os

args = sys.argv
arg_len = len(args)

if arg_len <2:
	print("\n**** Written by DJP, 24/10/17 in Python 3.4 ****\n")
	print("This script takes B2G annotation exported in 'WEGO' format (tab-delim by gene, with a header), and puts it into a from useable by TopGO")
	
	print("\n**** USAGE ****\n")
	print("python B2G_to_topGO.py [annot file from B2G]\n\n\n")
else:

	input_file_name = args[1]
	count_N = 0
	in_file = open(input_file_name)
	
	out_file = open(input_file_name + "_fortopgo.txt", "w")	

	for line in in_file:
		count_N = count_N + 1
		line = line.rstrip("\n")
		if count_N > 1:
			line = line.split("\t")
			gene_name = line[0]
			line_sp_c = 0
			line_stick = ""
			for el in line:
				line_sp_c = line_sp_c + 1
				if line_sp_c > 1:
					line_stick = line_stick + ", " + el
			
			line_stick = line_stick.lstrip(", ")
			line_stick = gene_name + "\t" + line_stick
			
			out_file.write(line_stick + "\n")
			#print(line)
			#print(line_stick)
			
	print("\nFinished, Maya.\n\n")