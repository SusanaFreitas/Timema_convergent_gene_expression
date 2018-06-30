#### readcount_sex_asex_randomiser.py

import sys
import os
import getopt
import numpy

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:N:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_file_name			  = "NOTHINGSET"
output_filename_base      = "NOTHINGSET"
N_times                   = "NOTHINGSET"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** readcount_sex_asex_randomiser.py | written by DJP 22/06/18, Lausanne to Edinburgh ****\n")
		print("\nThis script randomises readcount allocation between asex and sex columns within pair for each gene. Note the switch is kept the same switch for all tissues")
		print("\nExample:\n")
		print("Gene_name\tTbiSFRep1\tTbiSFRep2\tTbiSFRep3\tTteAFRep1\tTteAFRep2\tTteAFRep3")
		print("OG-1000\t\t497\t\t663\t\t936\t\t863\t\t1100\t\t1019")
		print("OG-1001\t\t835\t\t568\t\t703\t\t860\t\t706\t\t918")
		print("OG-1003\t\t264\t\t309\t\t409\t\t706\t\t719\t\t888")
		print("\nbecomes:\n")
		print("Gene_name\tTbiSFRep1\tTbiSFRep2\tTbiSFRep3\tTteAFRep1\tTteAFRep2\tTteAFRep3")
		print("OG-1000\t\t863\t\t1100\t\t1019\t\t497\t\t663\t\t936 \t\t## switched")
		print("OG-1001\t\t835\t\t568\t\t703\t\t860\t\t706\t\t918 \t\t## not switched")
		print("OG-1003\t\t706\t\t719\t\t888\t\t264\t\t309\t\t409 \t\t## switched")
		
		print("\nONLY USE 10sp_orth_readcounts.csv with this script, or I doubt it will give useful output...\n\n")

		print("\n**** USAGE **** \n")

		print("python readcount_sex_asex_randomiser.py -i 10sp_orth_readcounts.csv -N [Number of permuted files] -o [output base name]\n\n")
		 

		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-o'):
		output_filename_base = arg
	elif opt in ('-N'):
		N_times  = arg
	else:
		print("i dont know")
		sys.exit(2)

try:
	int(N_times)
	N_times = int(N_times)
except:
	print("\n\nValue for N_times is not an integer, exiting\n\n")
	sys.exit(2)


#### make output_dict

### mk output dir
output_dir_name = output_filename_base + "_ASRAND"
try:
	os.mkdir(output_dir_name)
except OSError as exc:
	print("\n********* WARNING ************\n" + output_dir_name + " already exists. files may be overwritten.")


#### read data into dicts


#Gene_name,
# 1,2,3    Tbi_SF_WB_Re1,Tbi_SF_WB_Re2,Tbi_SF_WB_Re3,
# 4,5,6    Tce_SF_WB_Re1,Tce_SF_WB_Re2,Tce_SF_WB_Re3,
# 7,8,9    Tcm_SF_WB_Re1,Tcm_SF_WB_Re2,Tcm_SF_WB_Re3,
# 10,11,12 Tdi_AF_WB_Re1,Tdi_AF_WB_Re2,Tdi_AF_WB_Re3,
# 13,14,15 Tge_AF_WB_Re1,Tge_AF_WB_Re2,Tge_AF_WB_Re3,
# 16,17,18 Tms_AF_WB_Re1,Tms_AF_WB_Re2,Tms_AF_WB_Re3,
# 19,20,21 Tpa_SF_WB_Re1,Tpa_SF_WB_Re2,Tpa_SF_WB_Re3,
# 22,23,24 Tps_SF_WB_Re1,Tps_SF_WB_Re2,Tps_SF_WB_Re3,
# 25,26,27 Tsi_AF_WB_Re1,Tsi_AF_WB_Re2,Tsi_AF_WB_Re3,
# 28,29,30 Tte_AF_WB_Re1,Tte_AF_WB_Re2,Tte_AF_WB_Re3,

# 31,32,33 Tbi_SF_RT_Re1,Tbi_SF_RT_Re2,Tbi_SF_RT_Re3,
# 34,35,36 Tbi_SF_LG_Re1,Tbi_SF_LG_Re2,Tbi_SF_LG_Re3,

# 37,38,39 Tce_SF_RT_Re1,Tce_SF_RT_Re2,Tce_SF_RT_Re3,
# 40,41,42 Tce_SF_LG_Re1,Tce_SF_LG_Re2,Tce_SF_LG_Re3,

# 43,44,45 Tcm_SF_RT_Re1,Tcm_SF_RT_Re2,Tcm_SF_RT_Re3,
# 46,47,48 Tcm_SF_LG_Re1,Tcm_SF_LG_Re2,Tcm_SF_LG_Re3,

# 49,50,51 Tdi_AF_RT_Re1,Tdi_AF_RT_Re2,Tdi_AF_RT_Re3,
# 52,53,54 Tdi_AF_LG_Re1,Tdi_AF_LG_Re2,Tdi_AF_LG_Re3,

# 55,56,57 Tge_AF_RT_Re1,Tge_AF_RT_Re2,Tge_AF_RT_Re3,
# 58,59,60 Tge_AF_LG_Re1,Tge_AF_LG_Re2,Tge_AF_LG_Re3,

# 61,62,63 Tms_AF_RT_Re1,Tms_AF_RT_Re2,Tms_AF_RT_Re3,
# 64,65,66 Tms_AF_LG_Re1,Tms_AF_LG_Re2,Tms_AF_LG_Re3,

# 67,68,69 Tpa_SF_RT_Re1,Tpa_SF_RT_Re2,Tpa_SF_RT_Re3,
# 70,71,72 Tpa_SF_LG_Re1,Tpa_SF_LG_Re2,Tpa_SF_LG_Re3,

# 73,74,75 Tps_SF_RT_Re1,Tps_SF_RT_Re2,Tps_SF_RT_Re3,
# 76,77,78 Tps_SF_LG_Re1,Tps_SF_LG_Re2,Tps_SF_LG_Re3,

# 79,80,81 Tsi_AF_RT_Re1,Tsi_AF_RT_Re2,Tsi_AF_RT_Re3,
# 82,83,84 Tsi_AF_LG_Re1,Tsi_AF_LG_Re2,Tsi_AF_LG_Re3,

# 85,86,87 Tte_AF_RT_Re1,Tte_AF_RT_Re2,Tte_AF_RT_Re3,
# 88,89,90 Tte_AF_LG_Re1,Tte_AF_LG_Re2,Tte_AF_LG_Re3





in_file = open(in_file_name)

gene_name_l = []

WB_Tbi_dict = {}
WB_Tce_dict = {}
WB_Tcm_dict = {}
WB_Tpa_dict = {}
WB_Tps_dict = {}

RT_Tbi_dict = {}
RT_Tce_dict = {}
RT_Tcm_dict = {}
RT_Tpa_dict = {}
RT_Tps_dict = {}

LG_Tbi_dict = {}
LG_Tce_dict = {}
LG_Tcm_dict = {}
LG_Tpa_dict = {}
LG_Tps_dict = {}


line_N = 0
for line in in_file:
	line_N = line_N + 1
	line = line.rstrip("")
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		gene_name = line[0]
		gene_name_l.append(gene_name)
			
		WB_Tbi_dict[gene_name] = [line[1] + ","  + line[2] + ","  + line[3],  line[28] + "," + line[29] + "," + line[30]]
		WB_Tce_dict[gene_name] = [line[4] + ","  + line[5] + ","  + line[6],  line[16] + "," + line[17] + "," + line[18]]
		WB_Tcm_dict[gene_name] = [line[7] + ","  + line[8] + ","  + line[9],  line[25] + "," + line[26] + "," + line[27]]
		WB_Tpa_dict[gene_name] = [line[19] + "," + line[20] + "," + line[21], line[13] + "," + line[14] + "," + line[15]]
		WB_Tps_dict[gene_name] = [line[22] + "," + line[23] + "," + line[24], line[10] + "," + line[11] + "," + line[12]]

		RT_Tbi_dict[gene_name] = [line[31] + "," + line[32] + "," + line[33], line[85] + "," + line[86] + "," + line[87]]
		RT_Tce_dict[gene_name] = [line[37] + "," + line[38] + "," + line[39], line[61] + "," + line[62] + "," + line[63]]
		RT_Tcm_dict[gene_name] = [line[43] + "," + line[44] + "," + line[45], line[79] + "," + line[80] + "," + line[81]]
		RT_Tpa_dict[gene_name] = [line[67] + "," + line[68] + "," + line[69], line[55] + "," + line[56] + "," + line[57]]
		RT_Tps_dict[gene_name] = [line[73] + "," + line[74] + "," + line[75], line[49] + "," + line[50] + "," + line[51]]

		LG_Tbi_dict[gene_name] = [line[34] + "," + line[35] + "," + line[37], line[88] + "," + line[89] + "," + line[90]]
		LG_Tce_dict[gene_name] = [line[40] + "," + line[41] + "," + line[42], line[64] + "," + line[65] + "," + line[66]]
		LG_Tcm_dict[gene_name] = [line[46] + "," + line[47] + "," + line[48], line[82] + "," + line[83] + "," + line[83]]
		LG_Tpa_dict[gene_name] = [line[70] + "," + line[71] + "," + line[72], line[58] + "," + line[59] + "," + line[60]]
		LG_Tps_dict[gene_name] = [line[76] + "," + line[77] + "," + line[78], line[52] + "," + line[53] + "," + line[54]]

in_file.close()

################################################################################################
##

out_header = "Gene_name," + \
	"Tbi_SF_WB_Re1,Tbi_SF_WB_Re2,Tbi_SF_WB_Re3,Tte_AF_WB_Re1,Tte_AF_WB_Re2,Tte_AF_WB_Re3," + \
	"Tce_SF_WB_Re1,Tce_SF_WB_Re2,Tce_SF_WB_Re3,Tms_AF_WB_Re1,Tms_AF_WB_Re2,Tms_AF_WB_Re3," + \
	"Tcm_SF_WB_Re1,Tcm_SF_WB_Re2,Tcm_SF_WB_Re3,Tsi_AF_WB_Re1,Tsi_AF_WB_Re2,Tsi_AF_WB_Re3," + \
	"Tpa_SF_WB_Re1,Tpa_SF_WB_Re2,Tpa_SF_WB_Re3,Tge_AF_WB_Re1,Tge_AF_WB_Re2,Tge_AF_WB_Re3," + \
	"Tps_SF_WB_Re1,Tps_SF_WB_Re2,Tps_SF_WB_Re3,Tdi_AF_WB_Re1,Tdi_AF_WB_Re2,Tdi_AF_WB_Re3," + \
	"Tbi_SF_RT_Re1,Tbi_SF_RT_Re2,Tbi_SF_RT_Re3,Tte_AF_RT_Re1,Tte_AF_RT_Re2,Tte_AF_RT_Re3," + \
	"Tce_SF_RT_Re1,Tce_SF_RT_Re2,Tce_SF_RT_Re3,Tms_AF_RT_Re1,Tms_AF_RT_Re2,Tms_AF_RT_Re3," + \
	"Tcm_SF_RT_Re1,Tcm_SF_RT_Re2,Tcm_SF_RT_Re3,Tsi_AF_RT_Re1,Tsi_AF_RT_Re2,Tsi_AF_RT_Re3," + \
	"Tpa_SF_RT_Re1,Tpa_SF_RT_Re2,Tpa_SF_RT_Re3,Tge_AF_RT_Re1,Tge_AF_RT_Re2,Tge_AF_RT_Re3," + \
	"Tps_SF_RT_Re1,Tps_SF_RT_Re2,Tps_SF_RT_Re3,Tdi_AF_RT_Re1,Tdi_AF_RT_Re2,Tdi_AF_RT_Re3," + \
	"Tbi_SF_LG_Re1,Tbi_SF_LG_Re2,Tbi_SF_LG_Re3,Tte_AF_LG_Re1,Tte_AF_LG_Re2,Tte_AF_LG_Re3," + \
	"Tce_SF_LG_Re1,Tce_SF_LG_Re2,Tce_SF_LG_Re3,Tms_AF_LG_Re1,Tms_AF_LG_Re2,Tms_AF_LG_Re3," + \
	"Tcm_SF_LG_Re1,Tcm_SF_LG_Re2,Tcm_SF_LG_Re3,Tsi_AF_LG_Re1,Tsi_AF_LG_Re2,Tsi_AF_LG_Re3," + \
	"Tpa_SF_LG_Re1,Tpa_SF_LG_Re2,Tpa_SF_LG_Re3,Tge_AF_LG_Re1,Tge_AF_LG_Re2,Tge_AF_LG_Re3," + \
	"Tps_SF_LG_Re1,Tps_SF_LG_Re2,Tps_SF_LG_Re3,Tdi_AF_LG_Re1,Tdi_AF_LG_Re2,Tdi_AF_LG_Re3"

#print(out_header)

################################################################################################
#### flip sex - asex order within pair for each genes. KEEP same switch for all tissues

switch_list = [0,1]

run_times = 0

out_file_N_flipped = open(output_filename_base + "_N_flipped.csv", "w")
out_file_N_flipped.write("file,Flipped_0_times,Flipped_1_times,Flipped_2_times,Flipped_3_times,Flipped_4_times,Flipped_5_times\n")


for n in range(0, N_times):
	
	
	output_filename = os.path.join(output_dir_name,output_filename_base + "_SArand_N_" + str(n + 1) + ".csv")
	output_file = open(output_filename, "w")
	
	output_file.write(out_header + "\n")
	run_times = run_times + 1
	
	N_flip_0 = 0
	N_flip_1 = 0
	N_flip_2 = 0
	N_flip_3 = 0
	N_flip_4 = 0
	N_flip_5 = 0
	
	for gene in gene_name_l:
		Tbi_order = numpy.random.choice(switch_list, 2, replace = False)
		Tce_order = numpy.random.choice(switch_list, 2, replace = False)
		Tcm_order = numpy.random.choice(switch_list, 2, replace = False)
		Tpa_order = numpy.random.choice(switch_list, 2, replace = False)
		Tps_order = numpy.random.choice(switch_list, 2, replace = False)
		
		N_flipped = Tbi_order[0] + Tce_order[0] + Tcm_order[0] + Tpa_order[0] + Tps_order[0]
		
		if N_flipped == 0:
			N_flip_0 = N_flip_0 + 1 
		if N_flipped == 1:
			N_flip_1 = N_flip_1 + 1 		
		if N_flipped == 2:
			N_flip_2 = N_flip_2 + 1 		
		if N_flipped == 3:
			N_flip_3 = N_flip_3 + 1 		
		if N_flipped == 4:
			N_flip_4 = N_flip_4 + 1 		
		if N_flipped == 5:
			N_flip_5 = N_flip_5 + 1 		
		
		out_list = []
		
		Tbi_WB_rec = WB_Tbi_dict.get(gene)
		out_list.append(Tbi_WB_rec[Tbi_order[0]])
		out_list.append(Tbi_WB_rec[Tbi_order[1]])
	
		Tce_WB_rec = WB_Tce_dict.get(gene)
		out_list.append(Tce_WB_rec[Tce_order[0]])
		out_list.append(Tce_WB_rec[Tce_order[1]])
	
		Tcm_WB_rec = WB_Tcm_dict.get(gene)
		out_list.append(Tcm_WB_rec[Tcm_order[0]])
		out_list.append(Tcm_WB_rec[Tcm_order[1]])
		
		Tpa_WB_rec = WB_Tpa_dict.get(gene)
		out_list.append(Tpa_WB_rec[Tpa_order[0]])
		out_list.append(Tpa_WB_rec[Tps_order[1]])
		
		Tps_WB_rec = WB_Tps_dict.get(gene)
		out_list.append(Tps_WB_rec[Tps_order[0]])
		out_list.append(Tps_WB_rec[Tps_order[1]])
		
		Tbi_RT_rec = RT_Tbi_dict.get(gene)
		out_list.append(Tbi_RT_rec[Tbi_order[0]])
		out_list.append(Tbi_RT_rec[Tbi_order[1]])
	
		Tce_RT_rec = RT_Tce_dict.get(gene)
		out_list.append(Tce_RT_rec[Tce_order[0]])
		out_list.append(Tce_RT_rec[Tce_order[1]])
	
		Tcm_RT_rec = RT_Tcm_dict.get(gene)
		out_list.append(Tcm_RT_rec[Tcm_order[0]])
		out_list.append(Tcm_RT_rec[Tcm_order[1]])
		
		Tpa_RT_rec = RT_Tpa_dict.get(gene)
		out_list.append(Tpa_RT_rec[Tpa_order[0]])
		out_list.append(Tpa_RT_rec[Tpa_order[1]])
		
		Tps_RT_rec = RT_Tps_dict.get(gene)
		out_list.append(Tps_RT_rec[Tps_order[0]])
		out_list.append(Tps_RT_rec[Tps_order[1]])
		
		Tbi_LG_rec = LG_Tbi_dict.get(gene)
		out_list.append(Tbi_LG_rec[Tbi_order[0]])
		out_list.append(Tbi_LG_rec[Tbi_order[1]])
	
		Tce_LG_rec = LG_Tce_dict.get(gene)
		out_list.append(Tce_LG_rec[Tce_order[0]])
		out_list.append(Tce_LG_rec[Tce_order[1]])
	
		Tcm_LG_rec = LG_Tcm_dict.get(gene)
		out_list.append(Tcm_LG_rec[Tcm_order[0]])
		out_list.append(Tcm_LG_rec[Tcm_order[1]])
		
		Tpa_LG_rec = LG_Tpa_dict.get(gene)
		out_list.append(Tpa_LG_rec[Tpa_order[0]])
		out_list.append(Tpa_LG_rec[Tpa_order[1]])
		
		Tps_LG_rec = LG_Tps_dict.get(gene)
		out_list.append(Tps_LG_rec[Tps_order[0]])
		out_list.append(Tps_LG_rec[Tps_order[1]])
		
		out_line = ""
		for i in out_list:
			out_line = out_line + "," + i
		
		#print(out_line)
		
		
		output_file.write(gene + out_line + "\n")
	
	out_file_N_flipped.write(output_filename_base + "_SArand_N_" + str(n + 1) + ".csv," + str(N_flip_0) + "," + str(N_flip_1) + "," + str(N_flip_2) + "," + str(N_flip_3) + "," + str(N_flip_4) + "," + str(N_flip_5) + "\n")
	
	output_file.close()

out_file_N_flipped.close()


print("\nFinished Lukas. Run " + str(run_times) + " permutations\n")
print("Output files are in " + output_dir_name + "\n\n") 