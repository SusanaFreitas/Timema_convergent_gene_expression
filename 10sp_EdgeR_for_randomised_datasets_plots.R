

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### libraries

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
library(modeest)

print (sessionInfo())

# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.5

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
 # [1] modeest_2.1         raster_2.5-8        sp_1.2-5            pvclust_2.0-0       vegan_2.4-4         permute_0.9-4       RColorBrewer_1.1-2  pheatmap_1.0.8      gtable_0.2.0       
# [10] stringr_1.2.0       cowplot_0.8.0       lattice_0.20-35     ggplot2_2.2.1       gridExtra_2.3       VennDiagram_1.6.17  futile.logger_1.4.3

# loaded via a namespace (and not attached):
 # [1] Rcpp_0.12.13         cluster_2.0.6        magrittr_1.5         MASS_7.3-47          munsell_0.4.3        colorspace_1.3-2     rlang_0.1.6          plyr_1.8.4           tools_3.4.1         
# [10] parallel_3.4.1       nlme_3.1-131         mgcv_1.8-22          lambda.r_1.2         lazyeval_0.2.1       tibble_1.3.4         Matrix_1.2-11        futile.options_1.0.0 stringi_1.1.5       
# [19] compiler_3.4.1       scales_0.5.0      


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# Data

dat1 <- read.csv("rand_Nconvergentgenes_out_tidied.csv")
head(dat1)
str(dat1)

N_WB_converge = 57
N_RT_converge = 203
N_LG_converge = 206


### stats

#### modes

N_conv_modes = c(mlv(dat1$WB, , method = "mfv")$M, mlv(dat1$RT, , method = "mfv")$M, mlv(dat1$LG, , method = "mfv")$M)


N_conv_modes_df <- as.data.frame(cbind(
c("WB", "RT", "LG"),
N_conv_modes,
c(N_WB_converge, N_RT_converge, N_LG_converge)
))

colnames(N_conv_modes_df) <- c("Tissue","Mode_of_perm","Obs")
N_conv_modes_df

#### prob of observing a value at least as high as I do: 
WB_p <- length(subset(dat1, dat1$WB >= N_WB_converge)[,1]) / 10000
RT_p <- length(subset(dat1, dat1$RT >= N_RT_converge)[,1]) / 10000
LG_p <- length(subset(dat1, dat1$LG >= N_LG_converge)[,1]) / 10000


## as 0 is incorrect, < 1/10000 is all I can say
if(WB_p == 0){
	WB_p = 1/10000
} 
if(RT_p == 0){
	RT_p = 1/10000
} 
if(LG_p == 0){
	LG_p = 1/10000
} 

max(dat1$WB)
max(dat1$RT)
max(dat1$LG)


### plot

P_WB <- ggplot(dat1, aes(x=WB)) +
    theme_bw() +
    geom_histogram(color="darkblue", fill="blue", binwidth=0.5,alpha = 0.2) +
    geom_vline(xintercept = N_WB_converge,  linetype = "longdash",color ="red") +
    geom_hline(yintercept = 0) +
    xlab("Number of genes showing convergent expression") +
    ylab("Number of permutations") +
	ggtitle(paste("Whole-body | 10000 permutations |", "p = ", WB_p)) + theme(plot.title = element_text(hjust = 0.5))

	
	
P_RT <- ggplot(dat1, aes(x=RT)) +
    theme_bw() +
    geom_histogram(color="darkblue", fill="blue", binwidth=0.5,alpha = 0.2) +
    geom_vline(xintercept = N_RT_converge,  linetype = "longdash",color ="red") +
    geom_hline(yintercept = 0) +
    xlab("Number of genes showing convergent expression") +
    ylab("Number of permutations") +
	ggtitle(paste("Reproductive tract | 10000 permutations |", "p < ", RT_p)) + theme(plot.title = element_text(hjust = 0.5))

	
	
P_LG <- ggplot(dat1, aes(x=LG)) +
    theme_bw() +
    geom_histogram(color="darkblue", fill="blue", binwidth=0.5,alpha = 0.2) +
    geom_vline(xintercept = N_LG_converge,  linetype = "longdash",color ="red") +
    geom_hline(yintercept = 0) +
    xlab("Number of genes showing convergent expression") +
    ylab("Number of permutations") +
	ggtitle(paste("Legs | 10000 permutations |", "p < ", LG_p)) + theme(plot.title = element_text(hjust = 0.5))

	

#### output

png(filename = "10sp_Noverall_eff_genes_perm_v2.png", width = 6, height = 15, units = "in", bg = "white", res = 300)
plot_grid(P_WB, P_RT, P_LG, ncol = 1, nrow = 3)
dev.off()







