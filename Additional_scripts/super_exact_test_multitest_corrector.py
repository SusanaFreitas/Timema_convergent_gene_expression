### super_exact_test_multitest_corrector.py

import statsmodels.stats.multitest as smm
import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


input_file_name = "NOTHINGSET"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** super_exact_test_multitest_corrector.py | Written by DJP, 24/11/17 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Corrects output of superexacttest for multiple tests using BH FDR \n") 
		
		print("***** USAGE **** \n") 
		print("python super_exact_test_multitest_corrector.py -i [input_file_name]\n")
		
		print("***** OPTIONS **** \n") 
		print("-i\t[input_file_name]\tinput file\t output from superexacttest via something like Cons_of_sex_asex_shifts_10sp_glmLRT.R\n\n")
		

		sys.exit(2)
		
	elif opt in ('-i'):
		input_file_name = arg
	else:
		print("i dont know")
		sys.exit(2)

if input_file_name == "NOTHINGSET":
	print("\n\nNo input_file_name, exiting! \n")
	sys.exit(2)



#### read input file to get dict of adj pvals

adjusted_pvals_dict = {}

set_name_list = []
pval_list = []

count_N = 0
in_file = open(input_file_name)
for line in in_file:
	count_N = count_N + 1
	line = line.rstrip("\n").split(",")
	
	if count_N != 1:
		if line[5] != 'NA':
			set_name_list.append(line[0])
			pval_list.append(float(line[5]))
in_file.close()

##### adjust pvals
rej, pval_corr = smm.multipletests(pval_list, method="fdr_bh")[:2]
pval_corr_list = list(pval_corr)

for i in range(0,len(set_name_list)):
	set_n = set_name_list[i]
	pc = pval_corr_list[i]
	adjusted_pvals_dict[set_n] = pc
	


#### output with FDR

outfilename = input_file_name + "_withFDR.csv"
outfile = open(outfilename, "w")

count_N = 0
in_file = open(input_file_name)
for line in in_file:
	count_N = count_N + 1
	start_line =  line.rstrip("\n").split(",")
	genes = line.rstrip("\n").split(",")
	header = ""
	if count_N == 1:
		start_line = start_line[0:6]
		start_line.append('"FDR_BH"')
		del genes[0:6]
		out_list = start_line + genes
		
		for j in out_list :
			header = header + "," + j
		header = header.lstrip(",")
		outfile.write(header + "\n")
	
	else:
		del genes[0:6]
		
		set_name = start_line[0]
		adj_pval = adjusted_pvals_dict.get(set_name)
		if str(adj_pval) == "None":
			adj_pval = "NA"
		
		
		start_line = start_line[0:6]
		start_line.append(adj_pval)
		out_list = start_line + genes
		
		out_line = ""
		for el in out_list:
			out_line = str(out_line) + "," + str(el)
		
		out_line = out_line.lstrip(",")
		

		outfile.write(out_line + "\n")

print("\nFinished, Sax\n")





