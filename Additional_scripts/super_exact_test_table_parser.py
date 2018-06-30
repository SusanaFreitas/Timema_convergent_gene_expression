# super_exact_test_table_parser.py

import sys
import decimal

#print(sys.argv) ### prints args

args = sys.argv ## to access args
arg_len = len(args)
if arg_len <3 :
	print("\n**** Written by DJP, 18/09/17 in Python 3.4 ****\n")
	print("This program takes the output file from SuperExactTest (in R) AFTER FDR correction with super_exact_test_multitest_corrector.py).") 
	print("It tidies up the file, and produces a matrix for the two way intersections (Note 2 files one with FDR and one without)")
	print("\n**** USAGE **** \n")
	print("python super_exact_test_table_parser.py [name of input file] [set_names] \n")
	print("\n**** USAGE OPTIONS ****\n")
	print("set_names\tSuperExactTest labels the sets as Set1, Set2 ... SetN. Here need to add the real names IN ORDER as a comma delim list\n")
else:

	SET_Table_name = args[1]
	set_names = args[2]
	set_names_co = set_names
	
	set_names = set_names.rstrip(",").split(",")
	print("\nExpecting to find " + str(len(set_names)) + " sets in the file.\n")


	### first check the number of Sets in the file
	SET_Table = open(SET_Table_name)
	line_N = 0
	set_seen_names = set()
	for line in SET_Table:
		line = line.rstrip("\n")
		line_N  = line_N + 1
		if line_N > 1:
			line = line.split(",")
			set_name = line[0].strip('"').split(" ")[0]
			set_seen_names.add(set_name)
	
	SET_Table.close()	
	
	found_sets_N = len(set_seen_names)
	
	if found_sets_N == len(set_names):
		print("Found " + str(found_sets_N) + " sets in the file.\n")
		print("Sets 1 to " + str(found_sets_N) + " are labelled in the following order: " + set_names_co + "\n")
		
	else:
		print("\nERROR, found " + str(found_sets_N) + " sets in the file, Exiting!\n")
		sys.exit(2)

	
	switch_dict = {}
	for i in range(1,found_sets_N + 1):
		dict_k = "Set" + str(i)
		realName = set_names[i-1]
		switch_dict[dict_k] = realName
		
	##### Ok so now want to make a matrix for 'level 2' sets
	
	# put in dict
	two_way_dict = {}
	two_way_dict_F = {}
	SET_Table = open(SET_Table_name)
	line_N = 0
	for line in SET_Table:
		line = line.rstrip("\n")
		line_N  = line_N + 1
		if line_N > 1:
			line = line.split(",")
			
			set_name = line[0].strip('"').split(" ")
			set_name_1 = ""
			for el in set_name:
				set_name_1 = set_name_1 + el
			
			degree = int(line[1])
			obser = line[2]
			exp = line[3]
			pval = line[5]
			FDR = line[6]
			
			if degree == 2:
				#print(set_name_1)
				exp = decimal.Decimal(exp).to_integral_value()
				two_way_dict[set_name_1] = obser + " (" + str(exp) + ") " + pval
				two_way_dict_F[set_name_1] = obser + " (" + str(exp) + ") " + FDR
	
	SET_Table.close()		
	
	
	# arrange and export
	out_file_2way_name = SET_Table_name + "_twoway_matrix.csv"
	out_file_2way = open(out_file_2way_name, "w")
	
	out_file_2way_name_F = SET_Table_name + "_twoway_matrix_WITH_FDR.csv"
	out_file_2way_F = open(out_file_2way_name_F, "w")
	
	header_line = ""
	for el in set_names:
		header_line = header_line + "," + el
	
	out_file_2way.write(header_line + "\n")
	out_file_2way_F.write(header_line + "\n")
	
	for i in range(1,found_sets_N + 1):
		rec_toget = ""
		rec_toget_F = ""
		for j in range(1,found_sets_N + 1):
			rec = ""
			rec_F = ""
			if i == j:
				rec = ""
				rec_F = ""
			elif i > j:
				set_key = "Set" + str(j) + "&" + "Set" + str(i)
				rec = two_way_dict.get(set_key)
				rec_F = two_way_dict_F.get(set_key)
				
			else:
				rec = ""
				rec_F = ""
			#print(j)
			#print(rec)
			rec_toget = rec_toget + "," + str(rec)
			rec_toget_F = rec_toget_F + "," + str(rec_F)
			
		rec_toget = set_names[i-1] + rec_toget
		rec_toget_F = set_names[i-1] + rec_toget_F
		#print(rec_toget)
		out_file_2way.write(rec_toget + "\n")
		out_file_2way_F.write(rec_toget_F + "\n")
	out_file_2way.close()
	out_file_2way_F.close()
	
	
	#### tidy up main file and output

	out_file_allway_name = SET_Table_name + "_allway_table_WITH_FDR.csv"
	out_file_allway = open(out_file_allway_name , "w")
	SET_Table = open(SET_Table_name)
	line_N = 0
	for line in SET_Table:
		line = line.rstrip("\n")
		line_N  = line_N + 1
		if line_N == 1:
			out_file_allway.write(line + "\n")
		else:
			line = line.split(",")
			set_name = line[0].strip('"').split(" ")
			#print(set_name)
			new_set_name = ""
			for el in set_name:
				if el != "&":
					ent = switch_dict.get(el)
					new_set_name = new_set_name + " & " + ent
					
			new_set_name = new_set_name.lstrip(" & ")
			#print(new_set_name)
			line.pop(0)
			
			out_line = ""
			for ele in line:
				out_line = out_line + "," + ele
			
			out_line = new_set_name + out_line
			out_file_allway.write(out_line+ "\n")

	
	print("\nFinished, Nigral.\n\n")
		
			
		