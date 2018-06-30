#### Get_GO_term_parent_and_child_overlap_adjuster.py 

# http://purl.obolibrary.org/obo/go.obo

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'g:i:c:p:d:o:Th')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


GO_obo_file_path = "NOTHINGSET" ### download the .obo flat-file that describes the terms in the database from http://purl.obolibrary.org/obo/go.obo
in_files = "NOTHINGSET" 
col_for_p = "NOTHINGSET"
p_cut = "NOTHINGSET"
set_IDs = "NOTHINGSET"
use_test_data = "NO"
outfile_prefix = "NOTHINGSET"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** Get_GO_term_parent_and_child_overlap_adjuster.py | Written by DJP, 26/11/17 in Python 3.5 in Lausanne, Swiss ****\n")
		print("\nThis code is to deal with the problem that topographically close go terms from different enrichement analyses will not overlap in a venn diagram.")
		print("To deal with this the code takes output from TopGO and clusters significant GO terms into clusters that are 1 link away from any other significant term")
		print("The code then outputs cluster names by set. If there are numerous GO terms in a set from a cluster they are appended with a number (e.g. clust1__2, clust1__3) so the number of GO terms represented as pre-clustering\n")
		print("\n**** USAGE ****\n")
		print("Get_GO_term_parent_and_child_overlap_adjuster.py -g [path to obo GO file] [Options]")
		print("\n**** OPTIONS ****\n")
		print("-g\tpath to obo GO file\tThis script requires the .obo flat-file that describes GO terms realtionships. Download from from http://purl.obolibrary.org/obo/go.obo")		
		print("-i\tin files\t\tcomma-delim list of files from TopGO (one per set). E.g. file1.csv,file2.csv,file3.csv")
		print("-c\tcol_for_p\t\tcolumn number containing pvalues to be used")
		print("-p\tp_cut\t\t\tp-value threshold. Default = 0.05")
		print("-d\tset_IDs\t\t\tComma delimited list to use as set names. Default OFF (Use filenames for setnames)")
		print("-o\toutput prefix\t\tprefix for all output files. Default: GO_overlap_adj_")
		print("-T\tuse_test_data\t\tSpecify if want to run the code with the internal stored data. For expliantion of the test data see accompanying R script: Get_GO_term_parent_and_child_overlap_adjuster_test_data_expl.R")

		
		print("\n**** EXAMPLE USAGE (using internal test data)****\n")
		print("python Get_GO_term_parent_and_child_overlap_adjuster.py  -g /path/to/obo_GO_file/go.obo -T\n\n") 
		sys.exit(2)
		
	elif opt in ('-g'):
		GO_obo_file_path = arg
	elif opt in ('-i'):
		in_files = arg
	elif opt in ('-c'):
		col_for_p = arg
	elif opt in ('-p'):
		p_cut = arg
	elif opt in ('-d'):
		set_IDs = arg
	elif opt in ('-T'):
		use_test_data = "YES"
	elif opt in ('-o'):
		outfile_prefix = arg
	else:
		print("i dont know")
		sys.exit(2)

if GO_obo_file_path == "NOTHINGSET":
	print("NO GO_obo_file_path provided, exiting\n")
	sys.exit(2)

if outfile_prefix == "NOTHINGSET":
	outfile_prefix = "GO_overlap_adj_"
	print("\nNO outfile_prefix set, using default: " + outfile_prefix + "\n")


if use_test_data == "NO":
	
	if in_files == "NOTHINGSET":
	  print("no in_files set, exiting. If you want to use test data specify -T\n")
	  sys.exit(2)
	
	if col_for_p == "NOTHINGSET":
	  print("NO col_for_p set, exiting\n")
	  sys.exit(2)
	
	if p_cut == "NOTHINGSET":
		p_cut = 0.05
		print("\nNo pval cutoff selected, Defaulting to 0.05\n")
	
	if set_IDs == "NOTHINGSET":
		print("\nNo set_IDs set, Defaulting to filename\n")
		in_files_l = in_files.rstrip(",").split(",")
		set_IDs_l  = in_files_l
	else:
		in_files_l = in_files.rstrip(",").split(",")
		set_IDs_l    = set_IDs.rstrip(",").split(",")
		if len(set_IDs_l) != len(in_files_l):
			print("The number of file in the file list does not match the number of setIDs, exiting\n")
			sys.exit(2)
		else:
			print("File names: " + in_files)
			print("Set names: " + set_IDs)

else:
	print("\nUSING test data\n")	


def getTerm(stream):
	block = []
	for line in stream:
		if line.strip() == "[Term]" or line.strip() == "[Typedef]":
			break
		else:
			if line.strip() != "":
				block.append(line.strip())
	return block


def parseTagValue(term):
	data = {}
	for line in term:
		tag = line.split(': ',1)[0]
		value = line.split(': ',1)[1]
		if tag not in data:
			data[tag] = []

		data[tag].append(value)

	return data

oboFile = open(GO_obo_file_path,'r')

#declare a blank dictionary
#keys are the goids
terms = {}

#skip the file header lines
getTerm(oboFile)

#infinite loop to go through the obo file.
#Breaks when the term returned is empty, indicating end of file
while 1:
	#get the term using the two parsing functions
	term = parseTagValue(getTerm(oboFile))
	if len(term) != 0:
		termID = term['id'][0]

	#only add to the structure if the term has a is_a tag
	#the is_a value contain GOID and term definition
	#we only want the GOID
		if 'is_a' in term:
			termParents = [p.split()[0] for p in term['is_a']]

			if termID not in terms:
				#each goid will have two arrays of parents and children
				terms[termID] = {'p':[],'c':[]}
	
			#append parents of the current term
			terms[termID]['p'] = termParents
	
		#for every parent term, add this current term as children
		for termParent in termParents:
			if termParent not in terms:
			  terms[termParent] = {'p':[],'c':[]}
			terms[termParent]['c'].append(termID)
	else:
		break
		
# parent = terms[want_GO ]['p']
# children = terms[want_GO ]['c']


####### read in data

set_dict = {}

if use_test_data == "NO":
	
	file_Nu = -1
	for f in in_files_l:
		infile = open(f)
		file_Nu = file_Nu + 1
		line_N = 0
		sig_set = set()
		for line in infile:
			line_N = line_N + 1
			line = line.rstrip("\n").split("\t")
			if line_N == 1:
				print("column selction sig terms from: " + str(line[int(col_for_p) + -1]))
			else:
				pval_use = decimal.Decimal(line[int(col_for_p) + -1])
				if pval_use <= p_cut:
					GO_term = line[0]
					sig_set.add(GO_term)
					#print (line)
					
		set_dict[set_IDs_l[file_Nu]] = sig_set		
			
		infile.close()

### problem set (first term in list is related to all)

if use_test_data == "YES":
	
	set_dict["set1"] = set(["GO:0010243","GO:0043279","GO:0002181", "GO:0008039", "GO:0001522","GO:0042689"]) ## has 2 terms from caffine
	set_dict["set2"] = set(["GO:0031000","GO:0002181", "GO:0046692", "GO:0007394","GO:0046845"])
	set_dict["set3"] = set(["GO:0071313","GO:0002181", "GO:0090090", "GO:0032543","GO:0008045"])

all_set = set()
for s in set_dict:
	cur_s = set_dict.get(s)
	all_set = all_set | cur_s

### make term clusters (just with parent and child terms)

all_GOS = set()
for el in terms:
	all_GOS.add(el)

clust_dict = {}
N_obso = 0
term_N = 0
for el in all_set:
	term_N = term_N + 1
	if el in all_GOS: ### deal with when a term is not in the annotation (obsolete terms)
		parent = terms[el]['p']
		children = terms[el ]['c']
	else:
		N_obso = N_obso + 1
		parent = el
		children = el
		
	parent_group_set = set()
	child_group_set = set()
	for p in parent:
		if p in all_set:
			parent_group_set.add(p)
				
	for c in children:
		if c in all_set:    
			child_group_set.add(c)

	pco_set =  parent_group_set | child_group_set
	pco_set.add(el)
	
	clust_N = "cluster_" + str(term_N)
	clust_dict[clust_N] = pco_set
			
print(N_obso)

#### merge overlapping clusts

Keep_clust_dict = {}


already_got = set()
for el in clust_dict:
	
	curr_vals = clust_dict.get(el)
	was_in = "NO"
	dup_GO = ""
	
	for v in curr_vals:
		if v not in already_got:
			was_in = "NO"
			
		else:
			was_in = "YES"
			dup_GO = v
			break
	
	if was_in == "NO":
		already_got = already_got | curr_vals
		Keep_clust_dict[el] = clust_dict.get(el)
	else:
		for n in Keep_clust_dict:
			rec = Keep_clust_dict.get(n)
			if v in rec:
				rec = rec | curr_vals
				Keep_clust_dict[n] = rec
				already_got = already_got | curr_vals

	#print(was_in + "\t" + dup_GO)


### output cluster key table ## with N terms in clust + make GO to clust dict


clustkey_out_file_name = outfile_prefix + "clustkey_out" + ".txt"
clustkey_out_file = open(clustkey_out_file_name, "w")
GO_to_clust = {}

for el in Keep_clust_dict:
	rec = Keep_clust_dict.get(el)
	rec_len = len(rec)
	out_rec = ""
	for r in rec:
		out_rec = out_rec + "," + r
	
	out_rec = out_rec.lstrip(",")
	
	clustkey_out_file.write(el + "\t" + str(rec_len) + "\t" +  out_rec + "\n")
	for g in rec:
		GO_to_clust[g] = el
		#print (g + "\t" + str(el))

############ make output for venns

### use cluster names. if 2 go terms come from same cluster in a set then the name of the cluster is put twice (clusters will have __x on if peresnt more than once)

for s in set_dict:
	cur_s = set_dict.get(s)
	cur_out_file_name = outfile_prefix + s + ".txt"

	cur_out_file = open(cur_out_file_name, "w")
	cur_out_file.write(s + "\n")
	seen = set()
	for_N_list = []
	clust_list_out = []
	for g in cur_s:
		clust_n = GO_to_clust.get(g)
		clust_list_out.append(clust_n)

		
	for i in clust_list_out:
		out_Name = ""
		
		if i not in seen:
			seen.add(i)
			out_Name = i
		else:
			for_N_list.append(i)
			N_ist = for_N_list.count(i)
			i = i + "__" + str(N_ist)
			out_Name = i
			
		cur_out_file.write(out_Name + "\n") 
	
	#print(clust_n)
	
	cur_out_file.close()


print("\n\nFinished, Diabloceratops\n\n")
	








