
### using command line args
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### aborts if no command line args provided
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

curr_filename <- args[1]
outbdirprefixname   <- args[2]



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### libraries

library("edgeR")
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)
library(stringr)
library(gtable)
library(pheatmap)
library(RColorBrewer)
require(vegan)
library(pvclust)
library(raster)
library("SuperExactTest")

print (sessionInfo())

# # R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.3

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] SuperExactTest_0.99.4 raster_2.5-8          sp_1.2-5              pvclust_2.0-0        
 # [5] vegan_2.4-4           permute_0.9-4         RColorBrewer_1.1-2    stringr_1.2.0        
 # [9] VennDiagram_1.6.17    futile.logger_1.4.3   edgeR_3.18.1          limma_3.32.10        
# [13] gtable_0.2.0          lattice_0.20-35       gridExtra_2.3         cowplot_0.8.0        
# [17] pheatmap_1.0.8        ggplot2_2.2.1        

# loaded via a namespace (and not attached):
 # [1] Rcpp_0.12.13         cluster_2.0.6        magrittr_1.5         MASS_7.3-47         
 # [5] munsell_0.4.3        colorspace_1.3-2     rlang_0.1.6          plyr_1.8.4          
 # [9] tools_3.4.1          parallel_3.4.1       nlme_3.1-131         mgcv_1.8-22         
# [13] lambda.r_1.2         lazyeval_0.2.1       tibble_1.3.4         Matrix_1.2-11       
# [17] futile.options_1.0.0 stringi_1.1.5        compiler_3.4.1       scales_0.5.0        
# [21] locfit_1.5-9.1  


print(curr_filename)
print(outbdirprefixname )

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# Data

rawdata_10sp <- read.csv(curr_filename, check.names=FALSE, stringsAsFactors=FALSE)

##### output
length(strsplit(strsplit(curr_filename, ".csv")[[1]][1], "/")[[1]])

outfilebasename   <- strsplit(strsplit(curr_filename, ".csv")[[1]][1], "/")[[1]][length(strsplit(strsplit(curr_filename, ".csv")[[1]][1], "/")[[1]])]
#print(outfilebasename )

dir.create(paste(outbdirprefixname, "_Nconvergentgenes_out", sep = ""), showWarnings = FALSE)
dir.create(paste(outbdirprefixname, "_genetables_out", sep = ""), showWarnings = FALSE)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
## into DGE structure 

#### Whole-body samples
y_WB_10sp_UF <- DGEList(counts=rawdata_10sp[,c(
"Tbi_SF_WB_Re1","Tbi_SF_WB_Re2","Tbi_SF_WB_Re3","Tte_AF_WB_Re1","Tte_AF_WB_Re2","Tte_AF_WB_Re3",
"Tce_SF_WB_Re1","Tce_SF_WB_Re2","Tce_SF_WB_Re3","Tms_AF_WB_Re1","Tms_AF_WB_Re2","Tms_AF_WB_Re3",
"Tcm_SF_WB_Re1","Tcm_SF_WB_Re2","Tcm_SF_WB_Re3","Tsi_AF_WB_Re1","Tsi_AF_WB_Re2","Tsi_AF_WB_Re3",
"Tpa_SF_WB_Re1","Tpa_SF_WB_Re2","Tpa_SF_WB_Re3","Tge_AF_WB_Re1","Tge_AF_WB_Re2","Tge_AF_WB_Re3",
"Tps_SF_WB_Re1","Tps_SF_WB_Re2","Tps_SF_WB_Re3","Tdi_AF_WB_Re1","Tdi_AF_WB_Re2","Tdi_AF_WB_Re3"
)], genes=rawdata_10sp[,1:1])


#### Reproductive tract samples
y_RT_10sp_UF <- DGEList(counts=rawdata_10sp[,c(
"Tbi_SF_RT_Re1","Tbi_SF_RT_Re2","Tbi_SF_RT_Re3","Tte_AF_RT_Re1","Tte_AF_RT_Re2","Tte_AF_RT_Re3",
"Tce_SF_RT_Re1","Tce_SF_RT_Re2","Tce_SF_RT_Re3","Tms_AF_RT_Re1","Tms_AF_RT_Re2","Tms_AF_RT_Re3",
"Tcm_SF_RT_Re1","Tcm_SF_RT_Re2","Tcm_SF_RT_Re3","Tsi_AF_RT_Re1","Tsi_AF_RT_Re2","Tsi_AF_RT_Re3",
"Tpa_SF_RT_Re1","Tpa_SF_RT_Re2","Tpa_SF_RT_Re3","Tge_AF_RT_Re1","Tge_AF_RT_Re2","Tge_AF_RT_Re3",
"Tps_SF_RT_Re1","Tps_SF_RT_Re2","Tps_SF_RT_Re3","Tdi_AF_RT_Re1","Tdi_AF_RT_Re2","Tdi_AF_RT_Re3"
)], genes=rawdata_10sp[,1:1])

#### Leg samples
y_LG_10sp_UF <- DGEList(counts=rawdata_10sp[,c(
"Tbi_SF_LG_Re1","Tbi_SF_LG_Re2","Tbi_SF_LG_Re3","Tte_AF_LG_Re1","Tte_AF_LG_Re2","Tte_AF_LG_Re3",
"Tce_SF_LG_Re1","Tce_SF_LG_Re2","Tce_SF_LG_Re3","Tms_AF_LG_Re1","Tms_AF_LG_Re2","Tms_AF_LG_Re3",
"Tcm_SF_LG_Re1","Tcm_SF_LG_Re2","Tcm_SF_LG_Re3","Tsi_AF_LG_Re1","Tsi_AF_LG_Re2","Tsi_AF_LG_Re3",
"Tpa_SF_LG_Re1","Tpa_SF_LG_Re2","Tpa_SF_LG_Re3","Tge_AF_LG_Re1","Tge_AF_LG_Re2","Tge_AF_LG_Re3",
"Tps_SF_LG_Re1","Tps_SF_LG_Re2","Tps_SF_LG_Re3","Tdi_AF_LG_Re1","Tdi_AF_LG_Re2","Tdi_AF_LG_Re3"
)], genes=rawdata_10sp[,1:1])



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (code)

filt_and_norm_maj <- function(y,cpm_cut,cut_in_Nsams){
	
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples
	head(cpm(y)) 
	keep <- 
	  rowSums(cpm(y[,1:3])> cpm_cut) >= cut_in_Nsams & 
      rowSums(cpm(y[,4:6])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,7:9])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,10:12])> cpm_cut) >= cut_in_Nsams &     
      rowSums(cpm(y[,13:15])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,16:18])> cpm_cut) >= cut_in_Nsams &      
      rowSums(cpm(y[,19:21])> cpm_cut) >= cut_in_Nsams &    
      rowSums(cpm(y[,22:24])> cpm_cut) >= cut_in_Nsams &     
      rowSums(cpm(y[,25:27])> cpm_cut) >= cut_in_Nsams &
      rowSums(cpm(y[,28:30])> cpm_cut) >= cut_in_Nsams
   

	y <- y[keep,]

	y$samples$lib.size <- colSums(y$counts) # if filter recalc lib size
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	
	return(y)

}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### filt and norm     

y_WB_10sp <- filt_and_norm_maj(y_WB_10sp_UF,0.5,2)
y_RT_10sp <- filt_and_norm_maj(y_RT_10sp_UF,0.5,2)
y_LG_10sp <- filt_and_norm_maj(y_LG_10sp_UF,0.5,2)


#################################################################################################
### get DE_gene_table (code)

## returns full table and sig DE genes as a vector
get_DE_genes <- function(fita,FDRa){
	TT1 = topTags(fita, n =3000000000)
	TT2 = TT1$table
	temp_sig <- subset(TT2, TT2$FDR <= FDRa)	
	sig_genes = temp_sig$genes
	sig_logFC = temp_sig$logFC
	
	N_sig_genes <- length(sig_genes)
	#cat("Number of sig genes: ", N_sig_genes )
	
	r_list <- list("table" = TT2, "S_gene_list" = sig_genes, "S_logFC_list" = sig_logFC )
	return(r_list)
}

####### design matrix
sp_tiss = factor(c(
"Tbi","Tbi","Tbi","Tbi","Tbi","Tbi",
"Tce","Tce","Tce","Tce","Tce","Tce",
"Tcm","Tcm","Tcm","Tcm","Tcm","Tcm",
"Tpa","Tpa","Tpa","Tpa","Tpa","Tpa",
"Tps","Tps","Tps","Tps","Tps","Tps"
))

prop_tiss = factor(c(
"S","S","S","A","A","A",
"S","S","S","A","A","A",
"S","S","S","A","A","A",
"S","S","S","A","A","A",
"S","S","S","A","A","A"
))

### design matrixes
design_WB_10sp <- model.matrix(~sp_tiss * prop_tiss )
design_RT_10sp <- model.matrix(~sp_tiss * prop_tiss )
design_LG_10sp <- model.matrix(~sp_tiss * prop_tiss )


###### Est dispersion (samples)

y_WB_10sp <- estimateDisp(y_WB_10sp, design_WB_10sp )
y_RT_10sp <- estimateDisp(y_RT_10sp, design_RT_10sp)
y_LG_10sp <- estimateDisp(y_LG_10sp, design_LG_10sp)


#####################################################################################################################################################################
###### fit glm 

fit_4inter_WB_10sp <- glmFit(y_WB_10sp, design_WB_10sp,robust=TRUE)
fit_4inter_RT_10sp <- glmFit(y_RT_10sp, design_RT_10sp,robust=TRUE)
fit_4inter_LG_10sp <- glmFit(y_LG_10sp, design_LG_10sp,robust=TRUE)

#### overall interaction

fit_4inter_WB_10sp_c <- glmLRT(fit_4inter_WB_10sp,coef=7:10)
fit_4inter_RT_10sp_c <- glmLRT(fit_4inter_RT_10sp,coef=7:10)
fit_4inter_LG_10sp_c <- glmLRT(fit_4inter_LG_10sp,coef=7:10)

TTT_WB_10sp_S_A_INTER <- get_DE_genes(fit_4inter_WB_10sp_c, 0.05)
TTT_RT_10sp_S_A_INTER <- get_DE_genes(fit_4inter_RT_10sp_c, 0.05)
TTT_LG_10sp_S_A_INTER <- get_DE_genes(fit_4inter_LG_10sp_c, 0.05)

####### export full tables

setwd(paste(outbdirprefixname, "_genetables_out", sep = ""))

write.table(TTT_WB_10sp_S_A_INTER$table, paste(outfilebasename, "TTT_WB_10sp_S_A_INTER.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_10sp_S_A_INTER$table, paste(outfilebasename, "TTT_RT_10sp_S_A_INTER.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_10sp_S_A_INTER$table, paste(outfilebasename, "TTT_LG_10sp_S_A_INTER.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)


##### get overall effects using additive model

### design matrix WB_match_tiss
design_WB_10sp_a <- model.matrix(~sp_tiss + prop_tiss )
rownames(design_WB_10sp_a) <- colnames(y_WB_10sp)
design_WB_10sp_a

### design matrix RT
design_RT_10sp_a <- model.matrix(~sp_tiss + prop_tiss )
rownames(design_RT_10sp_a) <- colnames(y_RT_10sp)
design_RT_10sp_a

### design matrix LG
design_LG_10sp_a <- model.matrix(~sp_tiss + prop_tiss )
rownames(design_LG_10sp_a) <- colnames(y_LG_10sp)
design_LG_10sp_a

fit_4oveeff_WB_10sp <- glmFit(y_WB_10sp, design_WB_10sp_a,robust=TRUE)
fit_4oveeff_RT_10sp <- glmFit(y_RT_10sp, design_RT_10sp_a,robust=TRUE)
fit_4oveeff_LG_10sp <- glmFit(y_LG_10sp, design_LG_10sp_a,robust=TRUE)

fit_4propeff_WB_10sp <- glmLRT(fit_4oveeff_WB_10sp,coef=6)
fit_4propeff_RT_10sp <- glmLRT(fit_4oveeff_RT_10sp,coef=6)
fit_4propeff_LG_10sp <- glmLRT(fit_4oveeff_LG_10sp,coef=6)

TTT_WB_10sp_S_A_propeff <- get_DE_genes(fit_4propeff_WB_10sp, 0.05)
TTT_RT_10sp_S_A_propeff <- get_DE_genes(fit_4propeff_RT_10sp, 0.05)
TTT_LG_10sp_S_A_propeff <- get_DE_genes(fit_4propeff_LG_10sp, 0.05)



####### export full tables

write.table(TTT_WB_10sp_S_A_propeff$table, paste(outfilebasename, "TTT_WB_10sp_S_A_propeff.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_RT_10sp_S_A_propeff$table, paste(outfilebasename, "TTT_RT_10sp_S_A_propeff.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)
write.table(TTT_LG_10sp_S_A_propeff$table, paste(outfilebasename, "TTT_LG_10sp_S_A_propeff.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)



### genes showing sig OVERALL S vs A effect

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}


Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}


WB_10sp_S_A_overall_eff_genes <- TTT_WB_10sp_S_A_propeff$S_gene_list
RT_10sp_S_A_overall_eff_genes <- TTT_RT_10sp_S_A_propeff$S_gene_list
LG_10sp_S_A_overall_eff_genes <- TTT_LG_10sp_S_A_propeff$S_gene_list

### genes with sig interaction

WB_10sp_S_A_INTER_genes <- TTT_WB_10sp_S_A_INTER$S_gene_list
RT_10sp_S_A_INTER_genes <- TTT_RT_10sp_S_A_INTER$S_gene_list
LG_10sp_S_A_INTER_genes <- TTT_LG_10sp_S_A_INTER$S_gene_list


WB_10sp_inter_DE_venn_overall_effSA_list <- list(WB_ove_eff = WB_10sp_S_A_overall_eff_genes, WB_inter = WB_10sp_S_A_INTER_genes)
RT_10sp_inter_DE_venn_overall_effSA_list <- list(RT_ove_eff = RT_10sp_S_A_overall_eff_genes, RT_inter = RT_10sp_S_A_INTER_genes)
LG_10sp_inter_DE_venn_overall_effSA_list <- list(LG_ove_eff = LG_10sp_S_A_overall_eff_genes, LG_inter = LG_10sp_S_A_INTER_genes)

WB_ove_eff_genes <- Setdiff(WB_10sp_inter_DE_venn_overall_effSA_list["WB_ove_eff"], WB_10sp_inter_DE_venn_overall_effSA_list["WB_inter"])
RT_ove_eff_genes <- Setdiff(RT_10sp_inter_DE_venn_overall_effSA_list["RT_ove_eff"], RT_10sp_inter_DE_venn_overall_effSA_list["RT_inter"])
LG_ove_eff_genes <- Setdiff(LG_10sp_inter_DE_venn_overall_effSA_list["LG_ove_eff"], LG_10sp_inter_DE_venn_overall_effSA_list["LG_inter"])
# length(WB_ove_eff_genes) 
# length(RT_ove_eff_genes) 
# length(LG_ove_eff_genes)


N_cons_genes_df <- as.data.frame(cbind(length(WB_ove_eff_genes), length(RT_ove_eff_genes), length(LG_ove_eff_genes)))
colnames(N_cons_genes_df ) <- c("WB","RT","LG")


write.table(WB_ove_eff_genes, paste(outfilebasename, "WB_convergent_genes.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)
write.table(RT_ove_eff_genes, paste(outfilebasename, "RT_convergent_genes.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)
write.table(LG_ove_eff_genes, paste(outfilebasename, "LG_convergent_genes.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)

setwd("../")
setwd(paste(outbdirprefixname, "_Nconvergentgenes_out", sep = ""))

write.table(N_cons_genes_df, paste(outfilebasename, "N_convergent_genes.csv", sep = "__"), sep = ',', quote = FALSE, row.names = FALSE)
