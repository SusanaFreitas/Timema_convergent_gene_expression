### topGO

### using command line args
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### aborts if no command line args provided
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

sp_use <- args[1]
#sp_use = "Tdi" ### test
print(sp_use)


# install
# source("http://bioconductor.org/biocLite.R") 
# biocLite() 
# source("http://bioconductor.org/biocLite.R")   
# biocLite("topGO")
# biocLite("ALL")
# biocLite("affyLib")

library(topGO)
library(ALL)
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
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
 # [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods  
# [10] base     

# other attached packages:
 # [1] SuperExactTest_0.99.4 ggplot2_2.2.1         gridExtra_2.3         VennDiagram_1.6.17   
 # [5] futile.logger_1.4.3   ALL_1.18.0            topGO_2.28.0          SparseM_1.77         
 # [9] GO.db_3.4.1           AnnotationDbi_1.38.2  IRanges_2.10.5        S4Vectors_0.14.7     
# [13] Biobase_2.36.2        graph_1.54.0          BiocGenerics_0.22.1  

# loaded via a namespace (and not attached):
 # [1] Rcpp_0.12.13         munsell_0.4.3        bit_1.1-12           colorspace_1.3-2    
 # [5] lattice_0.20-35      rlang_0.1.6          plyr_1.8.4           blob_1.1.0          
 # [9] gtable_0.2.0         DBI_0.7              lambda.r_1.2         matrixStats_0.52.2  
# [13] lazyeval_0.2.1       bit64_0.9-7          digest_0.6.12        tibble_1.3.4        
# [17] futile.options_1.0.0 memoise_1.1.0        RSQLite_2.0          compiler_3.4.1      
# [21] scales_0.5.0         pkgconfig_2.0.1   

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
#output

dir.create("crosssp_topGO_out/")
dir.create(paste("crosssp_topGO_out/to",sp_use, sep = ""))

#### load annotation from blast2GO.

## nr_annotated
geneID2GO_Tbi_nr <- readMappings(file = "Data/For_GOterm_anaylsis/tbi_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tce_nr <- readMappings(file = "Data/For_GOterm_anaylsis/tce_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tcm_nr <- readMappings(file = "Data/For_GOterm_anaylsis/tcm_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tpa_nr <- readMappings(file = "Data/For_GOterm_anaylsis/tpa_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")
geneID2GO_Tps_nr <- readMappings(file = "Data/For_GOterm_anaylsis/tps_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_nr.annot_fortopgo.txt")

## Droso annotated
geneID2GO_Tbi_DROSO <- readMappings(file = "Data/For_GOterm_anaylsis/tbi_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tce_DROSO <- readMappings(file = "Data/For_GOterm_anaylsis/tce_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tcm_DROSO <- readMappings(file = "Data/For_GOterm_anaylsis/tcm_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tpa_DROSO <- readMappings(file = "Data/For_GOterm_anaylsis/tpa_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")
geneID2GO_Tps_DROSO <- readMappings(file = "Data/For_GOterm_anaylsis/tps_rbbh_end_plus_medtrim3_200_rrna_free_arthropoda_mixed_noblasthit_cont+extra_for_b2g_blastx_DROSO.annot_fortopgo.txt")


###############################################################################################################################################
#### read in tables with genename and sig value
### Ranked gene lists contrasts (all to to current mapped ref)

# WB

RGL_Tbi_WB_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tbi_lrt_WB_match_tiss__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tce_WB_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tce_lrt_WB_match_tiss__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tcm_WB_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tcm_lrt_WB_match_tiss__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tpa_WB_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tpa_lrt_WB_match_tiss__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tps_WB_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tps_lrt_WB_match_tiss__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))

# RT
RGL_Tbi_RT_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tbi_lrt_RT__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tce_RT_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tce_lrt_RT__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tcm_RT_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tcm_lrt_RT__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tpa_RT_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tpa_lrt_RT__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tps_RT_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tps_lrt_RT__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))

# LG
RGL_Tbi_LG_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tbi_lrt_LG__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tce_LG_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tce_lrt_LG__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tcm_LG_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tcm_lrt_LG__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tpa_LG_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tpa_lrt_LG__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))
RGL_Tps_LG_10sp_S_A_contrasts <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/Tps_lrt_LG__single_sp_to",sp_use,"RBBH_Arth_NB_Mix_contrasts_rankedlist.txt",sep = "")))

### Ranked gene lists consistant (to current mapped ref)

# WB
RGL_tospref_WB_10sp_S_A_overall_eff <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/",sp_use,"_to",sp_use,"lrt_WB_match_tiss_10sp_S_A_propeff_overall_eff_rankedlist_interatbottom.txt",sep = "")))
                                                              
# RT
RGL_tospref_RT_10sp_S_A_overall_eff <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/",sp_use,"_to",sp_use,"lrt_RT_10sp_S_A_propeff_overall_eff_rankedlist_interatbottom.txt",sep = "")))

# LG
RGL_tospref_LG_10sp_S_A_overall_eff <- as.list(read.table(paste("Data/For_GOterm_anaylsis/cross_sp/to",sp_use,"/",sp_use,"_to",sp_use,"lrt_LG_10sp_S_A_propeff_overall_eff_rankedlist_interatbottom.txt",sep = "")))



#### need to get the genes as a named numeric vector

RGL_Tbi_WB_10sp_S_A_contrasts_GL <- RGL_Tbi_WB_10sp_S_A_contrasts$V2
RGL_Tce_WB_10sp_S_A_contrasts_GL <- RGL_Tce_WB_10sp_S_A_contrasts$V2
RGL_Tcm_WB_10sp_S_A_contrasts_GL <- RGL_Tcm_WB_10sp_S_A_contrasts$V2
RGL_Tpa_WB_10sp_S_A_contrasts_GL <- RGL_Tpa_WB_10sp_S_A_contrasts$V2
RGL_Tps_WB_10sp_S_A_contrasts_GL <- RGL_Tps_WB_10sp_S_A_contrasts$V2
RGL_Tbi_RT_10sp_S_A_contrasts_GL <- RGL_Tbi_RT_10sp_S_A_contrasts$V2
RGL_Tce_RT_10sp_S_A_contrasts_GL <- RGL_Tce_RT_10sp_S_A_contrasts$V2
RGL_Tcm_RT_10sp_S_A_contrasts_GL <- RGL_Tcm_RT_10sp_S_A_contrasts$V2
RGL_Tpa_RT_10sp_S_A_contrasts_GL <- RGL_Tpa_RT_10sp_S_A_contrasts$V2
RGL_Tps_RT_10sp_S_A_contrasts_GL <- RGL_Tps_RT_10sp_S_A_contrasts$V2
RGL_Tbi_LG_10sp_S_A_contrasts_GL <- RGL_Tbi_LG_10sp_S_A_contrasts$V2
RGL_Tce_LG_10sp_S_A_contrasts_GL <- RGL_Tce_LG_10sp_S_A_contrasts$V2
RGL_Tcm_LG_10sp_S_A_contrasts_GL <- RGL_Tcm_LG_10sp_S_A_contrasts$V2
RGL_Tpa_LG_10sp_S_A_contrasts_GL <- RGL_Tpa_LG_10sp_S_A_contrasts$V2
RGL_Tps_LG_10sp_S_A_contrasts_GL <- RGL_Tps_LG_10sp_S_A_contrasts$V2

RGL_tospref_WB_10sp_S_A_overall_eff_GL <- RGL_tospref_WB_10sp_S_A_overall_eff$V2
RGL_tospref_RT_10sp_S_A_overall_eff_GL <- RGL_tospref_RT_10sp_S_A_overall_eff$V2
RGL_tospref_LG_10sp_S_A_overall_eff_GL <- RGL_tospref_LG_10sp_S_A_overall_eff$V2


names(RGL_Tbi_WB_10sp_S_A_contrasts_GL) <- RGL_Tbi_WB_10sp_S_A_contrasts$V1
names(RGL_Tce_WB_10sp_S_A_contrasts_GL) <- RGL_Tce_WB_10sp_S_A_contrasts$V1
names(RGL_Tcm_WB_10sp_S_A_contrasts_GL) <- RGL_Tcm_WB_10sp_S_A_contrasts$V1
names(RGL_Tpa_WB_10sp_S_A_contrasts_GL) <- RGL_Tpa_WB_10sp_S_A_contrasts$V1
names(RGL_Tps_WB_10sp_S_A_contrasts_GL) <- RGL_Tps_WB_10sp_S_A_contrasts$V1
names(RGL_Tbi_RT_10sp_S_A_contrasts_GL) <- RGL_Tbi_RT_10sp_S_A_contrasts$V1
names(RGL_Tce_RT_10sp_S_A_contrasts_GL) <- RGL_Tce_RT_10sp_S_A_contrasts$V1
names(RGL_Tcm_RT_10sp_S_A_contrasts_GL) <- RGL_Tcm_RT_10sp_S_A_contrasts$V1
names(RGL_Tpa_RT_10sp_S_A_contrasts_GL) <- RGL_Tpa_RT_10sp_S_A_contrasts$V1
names(RGL_Tps_RT_10sp_S_A_contrasts_GL) <- RGL_Tps_RT_10sp_S_A_contrasts$V1
names(RGL_Tbi_LG_10sp_S_A_contrasts_GL) <- RGL_Tbi_LG_10sp_S_A_contrasts$V1
names(RGL_Tce_LG_10sp_S_A_contrasts_GL) <- RGL_Tce_LG_10sp_S_A_contrasts$V1
names(RGL_Tcm_LG_10sp_S_A_contrasts_GL) <- RGL_Tcm_LG_10sp_S_A_contrasts$V1
names(RGL_Tpa_LG_10sp_S_A_contrasts_GL) <- RGL_Tpa_LG_10sp_S_A_contrasts$V1
names(RGL_Tps_LG_10sp_S_A_contrasts_GL) <- RGL_Tps_LG_10sp_S_A_contrasts$V1

names(RGL_tospref_WB_10sp_S_A_overall_eff_GL) <- RGL_tospref_WB_10sp_S_A_overall_eff$V1
names(RGL_tospref_RT_10sp_S_A_overall_eff_GL) <- RGL_tospref_RT_10sp_S_A_overall_eff$V1
names(RGL_tospref_LG_10sp_S_A_overall_eff_GL) <- RGL_tospref_LG_10sp_S_A_overall_eff$V1

### run enrichment

setwd(paste("crosssp_topGO_out/to",sp_use, sep = ""))

run_enrichment <- function(genelist, ref, sig_for_genes, sig_for_GO){

	
	topDiffGenes <- function(allScore) {return(allScore < sig_for_genes)}
	ref = eval(parse(text=paste(ref,sep='')))
	
	#### make GOdata object
	#### setting node size as 5 so at least 5 genes must be annot per GO terms 
	#### do enrichment tests
	
	GODATA_BP = new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 5)

	### get N GOs used

	GO_term_use_BP_list = GODATA_BP@graph@nodes
	N_GO_term_use_BP = length(GODATA_BP@graph@nodes)
	result_GSEA_BP     <- runTest(GODATA_BP, algorithm = "elim", statistic = "ks")

	### combined tables
	allRes1_BP <- GenTable(GODATA_BP, GSEA = result_GSEA_BP, ranksOf = "GSEA", topNodes = length(GODATA_BP@graph@nodes), numChar = 200)
	sig_GSEA_BP_GO     = subset(allRes1_BP, allRes1_BP$GSEA < sig_for_GO)$GO.ID

	
	## return everything!
	out_list = list("N_GO_term_use_BP" = N_GO_term_use_BP, 
	                "GO_term_use_BP_list" = GO_term_use_BP_list, 
	                "allRes1_BP" = allRes1_BP, 
	                "sig_GSEA_BP_GO" = sig_GSEA_BP_GO,
	                "GODATA_BP" = GODATA_BP) 
	return(out_list)

}

#### run the enrichment stuff (0.05)

# contrasts nr

sp_use_2 = ifelse(sp_use == "Tbi", "Tbi", 
           ifelse(sp_use == "Tte", "Tbi",  
           ifelse(sp_use == "Tce", "Tce",  
           ifelse(sp_use == "Tms", "Tce",  
           ifelse(sp_use == "Tcm", "Tcm",             
           ifelse(sp_use == "Tsi", "Tcm",  
           ifelse(sp_use == "Tpa", "Tpa",  
           ifelse(sp_use == "Tge", "Tpa",  
           ifelse(sp_use == "Tps", "Tps", 
           ifelse(sp_use == "Tdi", "Tps", "ERROR"))))))))))       
           
           

RGL_Tbi_WB_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tbi_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tce_WB_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tce_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tcm_WB_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tcm_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tpa_WB_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tpa_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tps_WB_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tps_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tbi_RT_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tbi_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tce_RT_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tce_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tcm_RT_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tcm_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tpa_RT_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tpa_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tps_RT_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tps_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tbi_LG_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tbi_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tce_LG_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tce_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tcm_LG_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tcm_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tpa_LG_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tpa_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_Tps_LG_10sp_S_A_contrasts_nr_enrich = run_enrichment(RGL_Tps_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)

# contrasts DROSO
RGL_Tbi_WB_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tbi_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tce_WB_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tce_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tcm_WB_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tcm_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tpa_WB_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tpa_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tps_WB_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tps_WB_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tbi_RT_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tbi_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tce_RT_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tce_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tcm_RT_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tcm_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tpa_RT_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tpa_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tps_RT_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tps_RT_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tbi_LG_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tbi_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tce_LG_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tce_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tcm_LG_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tcm_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tpa_LG_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tpa_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_Tps_LG_10sp_S_A_contrasts_DROSO_enrich = run_enrichment(RGL_Tps_LG_10sp_S_A_contrasts_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)


# overall_eff nr
RGL_tospref_WB_10sp_S_A_overall_eff_nr_enrich = run_enrichment(RGL_tospref_WB_10sp_S_A_overall_eff_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_tospref_RT_10sp_S_A_overall_eff_nr_enrich = run_enrichment(RGL_tospref_RT_10sp_S_A_overall_eff_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)
RGL_tospref_LG_10sp_S_A_overall_eff_nr_enrich = run_enrichment(RGL_tospref_LG_10sp_S_A_overall_eff_GL, paste("geneID2GO_",sp_use_2,"_nr",sep = ""), 0.05, 0.05)


# overall_eff DROSO
RGL_tospref_WB_10sp_S_A_overall_eff_DROSO_enrich = run_enrichment(RGL_tospref_WB_10sp_S_A_overall_eff_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_tospref_RT_10sp_S_A_overall_eff_DROSO_enrich = run_enrichment(RGL_tospref_RT_10sp_S_A_overall_eff_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)
RGL_tospref_LG_10sp_S_A_overall_eff_DROSO_enrich = run_enrichment(RGL_tospref_LG_10sp_S_A_overall_eff_GL, paste("geneID2GO_",sp_use_2,"_DROSO",sep = ""), 0.05, 0.05)



###################################################################################################################################
#### How well did the annotation do?
############################

species = c(paste("to",sp_use,sep = ""))

Dr_N_BP_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_BP)

Dr_N_MF_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_MF)

Dr_N_CC_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_CC)

Dr_N_ALL_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_ALL)

Dr_N_BP_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_BP)

Dr_N_MF_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_MF)

Dr_N_CC_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_CC)

Dr_N_ALL_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_ALL)


Dr_N_BP_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_BP)

Dr_N_MF_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_MF)

Dr_N_CC_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_CC)

Dr_N_ALL_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_DROSO_enrich$N_GO_term_use_ALL)

nr_N_BP_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_BP)

nr_N_MF_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_MF)

nr_N_CC_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_CC)

nr_N_ALL_WB <- c(
RGL_tospref_WB_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_ALL)

nr_N_BP_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_BP)

nr_N_MF_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_MF)

nr_N_CC_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_CC)

nr_N_ALL_RT <- c(
RGL_tospref_RT_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_ALL)


nr_N_BP_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_BP)

nr_N_MF_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_MF)

nr_N_CC_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_CC)

nr_N_ALL_LG <- c(
RGL_tospref_LG_10sp_S_A_overall_eff_nr_enrich$N_GO_term_use_ALL)

NGOTERMS_BPGO_DROSOvsnr <- as.data.frame(cbind(species,Dr_N_BP_WB,nr_N_BP_WB,Dr_N_BP_RT,nr_N_BP_RT,Dr_N_BP_LG,nr_N_BP_LG))

## output tables

write.table(NGOTERMS_BPGO_DROSOvsnr,  paste("to",sp_use,"_NGOTERMS_BPGO_DROSOvsnr.csv",sep = ""), sep = ',', quote = FALSE, row.names = FALSE)

################## N sig Go



N_sig_GOs <- function(basename){

	# GSEA
	
	Tbi_N_GO_sig_GSEA = length(eval(parse(text=paste('RGL_Tbi_',basename,'$','sig_GSEA_ALL_GO',sep=''))))
	Tce_N_GO_sig_GSEA = length(eval(parse(text=paste('RGL_Tce_',basename,'$','sig_GSEA_ALL_GO',sep=''))))
	Tcm_N_GO_sig_GSEA = length(eval(parse(text=paste('RGL_Tcm_',basename,'$','sig_GSEA_ALL_GO',sep=''))))
	Tpa_N_GO_sig_GSEA = length(eval(parse(text=paste('RGL_Tpa_',basename,'$','sig_GSEA_ALL_GO',sep=''))))
	Tps_N_GO_sig_GSEA = length(eval(parse(text=paste('RGL_Tps_',basename,'$','sig_GSEA_ALL_GO',sep=''))))
	
	GSEA_vect = c(Tbi_N_GO_sig_GSEA,Tce_N_GO_sig_GSEA,Tcm_N_GO_sig_GSEA,Tpa_N_GO_sig_GSEA,Tps_N_GO_sig_GSEA)
	
	return(list("GSEA_vect" = GSEA_vect))

	
}





###################################################################################################################################
#### Venns
############################


five_sp_DE_venn <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2)

	return(venny.plot)
}


#### diff sp samples


## constrasts DROSO GSEA_BP_GO
RGL_WB_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
paste("to",sp_use,", whole-body, GSEA")
)

RGL_RT_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
paste("to",sp_use,", rep tract, GSEA")
)

RGL_LG_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_10sp_S_A_contrasts_DROSO_enrich$sig_GSEA_BP_GO,
paste("to",sp_use,", LEGS, GSEA")
)


# grid.arrange(gTree(children=RGL_WB_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO),
# gTree(children=RGL_RT_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO),
# gTree(children=RGL_LG_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO))

venn_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO <- arrangeGrob(gTree(children=RGL_WB_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file=paste("to",sp_use,"_venn_S_A_contrasts_DROSO_enrich_venn_GSEA_GO.png", sep = ""), venn_10sp_S_A_contrasts_DROSO_enrich_venn_GSEA_BP_GO, width = 6, height = 18)

#
## constrasts nr GSEA_BP_GO
RGL_WB_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_WB_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_WB_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_WB_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_WB_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_WB_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
paste("to",sp_use,", whole-body, GSEA")
)

RGL_RT_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_RT_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_RT_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_RT_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_RT_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_RT_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
paste("to",sp_use,", rep tract, GSEA")
)

RGL_LG_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO <- five_sp_DE_venn(
RGL_Tbi_LG_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tce_LG_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tcm_LG_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tpa_LG_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
RGL_Tps_LG_10sp_S_A_contrasts_nr_enrich$sig_GSEA_BP_GO,
paste("to",sp_use,", LEGS, GSEA")
)


grid.arrange(gTree(children=RGL_WB_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_RT_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO),
gTree(children=RGL_LG_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO))

# venn_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO <- arrangeGrob(gTree(children=RGL_WB_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO),
# gTree(children=RGL_RT_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO),
# gTree(children=RGL_LG_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO), ncol=1)

ggsave(file=paste("to",sp_use,"_venn_S_A_contrasts_nr_enrich_venn_GSEA_GO.png", sep = ""), venn_10sp_S_A_contrasts_nr_enrich_venn_GSEA_BP_GO, width = 6, height = 18)


#### superexact test

run_superexact_GOs_BP <- function(basename){
	
	# GSEA
	GSEA_Tbi_GO_sig = eval(parse(text=paste('RGL_Tbi_',basename,'$','sig_GSEA_BP_GO',sep='')))
	GSEA_Tce_GO_sig = eval(parse(text=paste('RGL_Tce_',basename,'$','sig_GSEA_BP_GO',sep='')))
	GSEA_Tcm_GO_sig = eval(parse(text=paste('RGL_Tcm_',basename,'$','sig_GSEA_BP_GO',sep='')))
	GSEA_Tpa_GO_sig = eval(parse(text=paste('RGL_Tpa_',basename,'$','sig_GSEA_BP_GO',sep='')))	
	GSEA_Tps_GO_sig = eval(parse(text=paste('RGL_Tps_',basename,'$','sig_GSEA_BP_GO',sep='')))

	GSEA_all_GO_sig <- list(GSEA_Tbi_GO_sig,GSEA_Tce_GO_sig,GSEA_Tcm_GO_sig,GSEA_Tpa_GO_sig,GSEA_Tps_GO_sig)

	GSEA_Tbi_GO_allused = eval(parse(text=paste('RGL_Tbi_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tce_GO_allused = eval(parse(text=paste('RGL_Tce_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tcm_GO_allused = eval(parse(text=paste('RGL_Tcm_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tpa_GO_allused = eval(parse(text=paste('RGL_Tpa_',basename,'$','GO_term_use_BP_list',sep='')))	
	GSEA_Tps_GO_allused = eval(parse(text=paste('RGL_Tps_',basename,'$','GO_term_use_BP_list',sep='')))	
	
	GSEA_all_GO_allused <- list(GSEA_Tbi_GO_allused,GSEA_Tce_GO_allused,GSEA_Tcm_GO_allused,GSEA_Tpa_GO_allused,GSEA_Tps_GO_allused)	

	## lens of GOs differ to make more conservative I will use the size of the GO intersection of BP GOs
	GSEA_interlength <- length(Intersect(GSEA_all_GO_allused))
	
	print(GSEA_interlength)
	
	sup_out_GSEA = summary(supertest(GSEA_all_GO_sig, n= GSEA_interlength))

	out_list = list("sup_out_GSEA" = sup_out_GSEA)
	return(out_list)

}

### run ## not corrected for multi tests

sup_WB_10sp_S_A_contrasts_DROSO_enrich_GSEA_BP <- run_superexact_GOs_BP("WB_10sp_S_A_contrasts_DROSO_enrich")$sup_out_GSEA
sup_RT_10sp_S_A_contrasts_DROSO_enrich_GSEA_BP <- run_superexact_GOs_BP("RT_10sp_S_A_contrasts_DROSO_enrich")$sup_out_GSEA
sup_LG_10sp_S_A_contrasts_DROSO_enrich_GSEA_BP <- run_superexact_GOs_BP("LG_10sp_S_A_contrasts_DROSO_enrich")$sup_out_GSEA

sup_WB_10sp_S_A_contrasts_nr_enrich_GSEA_BP <- run_superexact_GOs_BP("WB_10sp_S_A_contrasts_nr_enrich")$sup_out_GSEA
sup_RT_10sp_S_A_contrasts_nr_enrich_GSEA_BP <- run_superexact_GOs_BP("RT_10sp_S_A_contrasts_nr_enrich")$sup_out_GSEA
sup_LG_10sp_S_A_contrasts_nr_enrich_GSEA_BP <- run_superexact_GOs_BP("LG_10sp_S_A_contrasts_nr_enrich")$sup_out_GSEA

## output

write.csv(sup_WB_10sp_S_A_contrasts_DROSO_enrich_GSEA_BP$Table, file=paste("to",sp_use,"_sup_WB_S_A_contrasts_DROSO_enrich_GSEA_BP.csv", sep = ""), row.names=FALSE)
write.csv(sup_RT_10sp_S_A_contrasts_DROSO_enrich_GSEA_BP$Table, file=paste("to",sp_use,"_sup_RT_S_A_contrasts_DROSO_enrich_GSEA_BP.csv", sep = ""), row.names=FALSE)
write.csv(sup_LG_10sp_S_A_contrasts_DROSO_enrich_GSEA_BP$Table, file=paste("to",sp_use,"_sup_LG_S_A_contrasts_DROSO_enrich_GSEA_BP.csv", sep = ""), row.names=FALSE)

write.csv(sup_WB_10sp_S_A_contrasts_nr_enrich_GSEA_BP$Table, file=paste("to",sp_use,"_sup_WB_S_A_contrasts_nr_enrich_GSEA_BP.csv", sep = ""), row.names=FALSE)
write.csv(sup_RT_10sp_S_A_contrasts_nr_enrich_GSEA_BP$Table, file=paste("to",sp_use,"_sup_RT_S_A_contrasts_nr_enrich_GSEA_BP.csv", sep = ""), row.names=FALSE)
write.csv(sup_LG_10sp_S_A_contrasts_nr_enrich_GSEA_BP$Table, file=paste("to",sp_use,"_sup_LG_S_A_contrasts_nr_enrich_GSEA_BP.csv", sep = ""), row.names=FALSE)



#### GO enrichment tables - export

## overall

write.table(RGL_tospref_WB_10sp_S_A_overall_eff_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_WB_10sp_S_A_conv_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_tospref_RT_10sp_S_A_overall_eff_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_RT_10sp_S_A_conv_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_tospref_LG_10sp_S_A_overall_eff_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_LG_10sp_S_A_conv_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_tospref_WB_10sp_S_A_overall_eff_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_WB_10sp_S_A_conv_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_tospref_RT_10sp_S_A_overall_eff_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_RT_10sp_S_A_conv_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_tospref_LG_10sp_S_A_overall_eff_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_LG_10sp_S_A_conv_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)



## contrasts

write.table(RGL_Tbi_WB_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tbi_WB_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tce_WB_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tcm_WB_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tpa_WB_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tps_WB_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_RT_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tbi_RT_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tce_RT_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tcm_RT_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tpa_RT_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tps_RT_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_LG_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tbi_LG_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tce_LG_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tcm_LG_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tpa_LG_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_10sp_S_A_contrasts_DROSO_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tps_LG_10sp_S_A_contrasts_DROSO_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)


write.table(RGL_Tbi_WB_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tbi_WB_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_WB_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tce_WB_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_WB_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tcm_WB_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_WB_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tpa_WB_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_WB_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tps_WB_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_RT_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tbi_RT_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_RT_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tce_RT_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_RT_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tcm_RT_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_RT_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tpa_RT_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_RT_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tps_RT_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)

write.table(RGL_Tbi_LG_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tbi_LG_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tce_LG_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tce_LG_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tcm_LG_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tcm_LG_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tpa_LG_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tpa_LG_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Tps_LG_10sp_S_A_contrasts_nr_enrich$allRes1_BP, paste("to",sp_use,"_RGL_Tps_LG_10sp_S_A_contrasts_nr_enrich_allRes1_BP.txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE)

########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), paste("to",sp_use,"cross_sp_topGO.R_sessionInfo.txt", sep = ""))


