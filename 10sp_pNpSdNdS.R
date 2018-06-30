## 10sp_pNpSdNdS.R
### libs

library(ggplot2)
library(cowplot)
library(RColorBrewer)

print (sessionInfo())

# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.3

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] RColorBrewer_1.1-2 cowplot_0.8.0      ggplot2_2.2.1     

# loaded via a namespace (and not attached):
 # [1] colorspace_1.3-2 scales_0.5.0     compiler_3.4.1   lazyeval_0.2.1   plyr_1.8.4       gtable_0.2.0     tibble_1.3.4    
 # [8] Rcpp_0.12.13     grid_3.4.1       rlang_0.1.6      munsell_0.4.3   


### data

dat1 = read.csv("Data/dnds_pnps/10_sp_dnds_pnps.csv")
length(dat1[,1])
head(dat1)
dat1$cons_gene_by_tiss_ord = ordered(dat1$cons_gene_by_tiss, levels = c("YES", "NO") )

dir.create("10sp_dndspnps_output")

### split by tissue

dat1_WB = subset(dat1, dat1$tissue == "WB")
dat1_RT = subset(dat1, dat1$tissue == "RT")
dat1_LG = subset(dat1, dat1$tissue == "LG")


### split by species

dat1_WB_Tte = subset(dat1_WB, dat1_WB$sp_name == "Tte")
dat1_WB_Tms = subset(dat1_WB, dat1_WB$sp_name == "Tms")
dat1_WB_Tsi = subset(dat1_WB, dat1_WB$sp_name == "Tsi")
dat1_WB_Tge = subset(dat1_WB, dat1_WB$sp_name == "Tge")
dat1_WB_Tdi = subset(dat1_WB, dat1_WB$sp_name == "Tdi")

dat1_RT_Tte = subset(dat1_RT, dat1_RT$sp_name == "Tte")
dat1_RT_Tms = subset(dat1_RT, dat1_RT$sp_name == "Tms")
dat1_RT_Tsi = subset(dat1_RT, dat1_RT$sp_name == "Tsi")
dat1_RT_Tge = subset(dat1_RT, dat1_RT$sp_name == "Tge")
dat1_RT_Tdi = subset(dat1_RT, dat1_RT$sp_name == "Tdi")

dat1_LG_Tte = subset(dat1_LG, dat1_LG$sp_name == "Tte")
dat1_LG_Tms = subset(dat1_LG, dat1_LG$sp_name == "Tms")
dat1_LG_Tsi = subset(dat1_LG, dat1_LG$sp_name == "Tsi")
dat1_LG_Tge = subset(dat1_LG, dat1_LG$sp_name == "Tge")
dat1_LG_Tdi = subset(dat1_LG, dat1_LG$sp_name == "Tdi")


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### boxplots by tiss_by_species (I.e. classification is on each sp independently)


# dnds
plot_box_dnds = function(df){
	exp_labels <- c("Convergent", "Background" )

	P1 <- ggplot(df, aes(cons_gene_by_tiss_ord  , dnds)) + 
		#theme_bw() +
		geom_boxplot(aes(fill = factor(cons_gene_by_tiss_ord ))) + 
		theme(legend.position="none") +
		ggtitle(deparse(substitute(df))) +
		ylim(0,max(dat1$dnds, na.rm = TRUE) + 0.1)
	P2 = P1 + scale_fill_manual(values=c("#1b9e77","#7570b3")) + scale_x_discrete(labels= exp_labels) + 
		labs(x="Expression change group", y = "dN/dS")
	return(P2)
}


dnds_WB <- plot_box_dnds(dat1_WB)
dnds_RT <- plot_box_dnds(dat1_RT)
dnds_LG <- plot_box_dnds(dat1_LG)

dnds_WB_Tte <- plot_box_dnds(dat1_WB_Tte)
dnds_WB_Tms <- plot_box_dnds(dat1_WB_Tms)
dnds_WB_Tsi <- plot_box_dnds(dat1_WB_Tsi)
dnds_WB_Tge <- plot_box_dnds(dat1_WB_Tge)
dnds_WB_Tdi <- plot_box_dnds(dat1_WB_Tdi)

dnds_RT_Tte <- plot_box_dnds(dat1_RT_Tte)
dnds_RT_Tms <- plot_box_dnds(dat1_RT_Tms)
dnds_RT_Tsi <- plot_box_dnds(dat1_RT_Tsi)
dnds_RT_Tge <- plot_box_dnds(dat1_RT_Tge)
dnds_RT_Tdi <- plot_box_dnds(dat1_RT_Tdi)

dnds_LG_Tte <- plot_box_dnds(dat1_LG_Tte)
dnds_LG_Tms <- plot_box_dnds(dat1_LG_Tms)
dnds_LG_Tsi <- plot_box_dnds(dat1_LG_Tsi)
dnds_LG_Tge <- plot_box_dnds(dat1_LG_Tge)
dnds_LG_Tdi <- plot_box_dnds(dat1_LG_Tdi)


dnds_all <- plot_grid(
dnds_WB_Tte, dnds_WB_Tms, dnds_WB_Tdi, dnds_WB_Tsi, dnds_WB_Tge,
dnds_RT_Tte, dnds_RT_Tms, dnds_RT_Tdi, dnds_RT_Tsi, dnds_RT_Tge,
dnds_LG_Tte, dnds_LG_Tms, dnds_LG_Tdi, dnds_LG_Tsi, dnds_LG_Tge 
, ncol = 5, nrow = 3)

dnds_title <- ggdraw() + draw_label("dN/dS", fontface='bold')

png(filename = "10sp_dndspnps_output/dnds_all.png", width = 13, height = 13, units = "in", bg = "white", res = 300)
plot_grid(dnds_title, dnds_all, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
dev.off()


# pnps

plot_box_pnps = function(df){
	exp_labels <- c("Convergent", "None" )

	P1 <- ggplot(df, aes(cons_gene_by_tiss_ord, pnps)) + 
		#theme_bw() +
		geom_boxplot(aes(fill = factor(cons_gene_by_tiss_ord))) + 
		theme(legend.position="none") +
		ggtitle(deparse(substitute(df))) + 
		ylim(0,max(dat1$pnps, na.rm = TRUE) + 0.1)
	P2 = P1 + scale_fill_manual(values=c("#1b9e77","#7570b3")) + scale_x_discrete(labels= exp_labels) + 
		labs(x="Expression change group", y = "pN/pS")
	return(P2)
}


pnps_WB <- plot_box_pnps(dat1_WB)
pnps_RT <- plot_box_pnps(dat1_RT)
pnps_LG <- plot_box_pnps(dat1_LG)


pnps_WB_Tte <- plot_box_pnps(dat1_WB_Tte)
pnps_WB_Tms <- plot_box_pnps(dat1_WB_Tms)
pnps_WB_Tsi <- plot_box_pnps(dat1_WB_Tsi)
pnps_WB_Tge <- plot_box_pnps(dat1_WB_Tge)
pnps_WB_Tdi <- plot_box_pnps(dat1_WB_Tdi)

pnps_RT_Tte <- plot_box_pnps(dat1_RT_Tte)
pnps_RT_Tms <- plot_box_pnps(dat1_RT_Tms)
pnps_RT_Tsi <- plot_box_pnps(dat1_RT_Tsi)
pnps_RT_Tge <- plot_box_pnps(dat1_RT_Tge)
pnps_RT_Tdi <- plot_box_pnps(dat1_RT_Tdi)

pnps_LG_Tte <- plot_box_pnps(dat1_LG_Tte)
pnps_LG_Tms <- plot_box_pnps(dat1_LG_Tms)
pnps_LG_Tsi <- plot_box_pnps(dat1_LG_Tsi)
pnps_LG_Tge <- plot_box_pnps(dat1_LG_Tge)
pnps_LG_Tdi <- plot_box_pnps(dat1_LG_Tdi)



pnps_all <- plot_grid(
pnps_WB_Tte, pnps_WB_Tms, pnps_WB_Tdi, pnps_WB_Tsi, pnps_WB_Tge, 
pnps_RT_Tte, pnps_RT_Tms, pnps_RT_Tdi, pnps_RT_Tsi, pnps_RT_Tge,
pnps_LG_Tte, pnps_LG_Tms, pnps_LG_Tdi, pnps_LG_Tsi, pnps_LG_Tge
, ncol = 5, nrow = 3)

pnps_title <- ggdraw() + draw_label("pN/pS", fontface='bold')
# plot_grid(pnps_title, pnps_all, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
png(filename = "10sp_dndspnps_output/pnps_all.png", width = 13, height = 13, units = "in", bg = "white", res = 300)
plot_grid(pnps_title, pnps_all, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
dev.off()
getwd() ## where has my plot 


####### overall plot

png(filename = "10sp_dndspnps_output/dnds_pnps_overall_groupbytissueandsp.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
dnds_WB,dnds_RT,dnds_LG,
pnps_WB,pnps_RT,pnps_LG
, ncol = 3, nrow = 2)
dev.off()
getwd() ## where has my plot 



###################### significant? DO permutation t-test

### randomise 

## randomise data by assiging convergent status to genes randomly

rand_gene_wil <- function(df){
	uniq_g_list = as.character(unique(df$OG_name))
	
	## how many converg genes?
	
	df_con = subset(df, df$cons_gene_by_tiss == "YES")
	con_uniq_gene = length(unique(df_con$OG_name))
	print(con_uniq_gene)
	
	### randomly select 'conserved' genes
	
	r_samp <- sample(uniq_g_list, con_uniq_gene)
	#print(r_samp)
	
	rand_cons_stat = c()
	for(i in 1:length(df$OG_name)){
		gene_name <- as.character(df[i,2])
		#print(i)
		#print(gene_name)
		c_stat <- ifelse(gene_name %in% r_samp, "YES", "NO")
		#print(c_stat)
		rand_cons_stat = c(rand_cons_stat, c_stat)		
		}
	
	df_out = as.data.frame(cbind(as.character(df$OG_name), rand_cons_stat, df$dnds, df$pnps))
	df_out$V3 = as.numeric(as.character(df_out$V3)) 
	df_out$V4 = as.numeric(as.character(df_out$V4))
	
	colnames(df_out) = c("OG_name", "rand_cons_stat", "dnds", "pnps") 
	return(df_out)
}

###### ttest

loop_t <- function(df, N_times){

	dnds_t_stat  = c()
	pnps_t_stat  = c()
	
	for(i in 1:N_times){
		print(i)	
		curr_df = rand_gene_wil(df)
		curr_df_cons = subset(curr_df, curr_df$rand_cons_stat == "YES")
		curr_df_Notcons = subset(curr_df, curr_df$rand_cons_stat == "NO")


		dnds_tt <- t.test(curr_df_cons$dnds, curr_df_Notcons$dnds, paired = FALSE, alternative = c("two.sided"))		
		pnps_tt <- t.test(curr_df_cons$pnps, curr_df_Notcons$pnps, paired = FALSE, alternative = c("two.sided"))				
		dnds_tt_W <- as.vector(dnds_tt$statistic)
		pnps_tt_W <- as.vector(pnps_tt$statistic)
		dnds_t_stat = c(dnds_t_stat, dnds_tt_W)
		pnps_t_stat = c(pnps_t_stat, pnps_tt_W )

	}
	
	out_list = list("dnds_t_stats" = dnds_t_stat, "pnps_t_stats" = pnps_t_stat)	
	return(out_list)

}


WB_t_rand <- loop_t(dat1_WB, 10000)
WB_t_rand_df <- as.data.frame(cbind(WB_t_rand$dnds_t_stats, WB_t_rand$pnps_t_stats))
colnames(WB_t_rand_df) <- c("dnds_t_teststat","pnps_t_teststat")

RT_t_rand <- loop_t(dat1_RT, 10000)
RT_t_rand_df <- as.data.frame(cbind(RT_t_rand$dnds_t_stats, RT_t_rand$pnps_t_stats))
colnames(RT_t_rand_df) <- c("dnds_t_teststat","pnps_t_teststat")

LG_t_rand <- loop_t(dat1_LG, 10000)
LG_t_rand_df <- as.data.frame(cbind(LG_t_rand$dnds_t_stats, LG_t_rand$pnps_t_stats))
colnames(LG_t_rand_df) <- c("dnds_t_teststat","pnps_t_teststat")


#### get observed test statistics


ty_tests_2 = function(df){
	
	df_name = deparse(substitute(df))
	sp_use = strsplit(df_name, "_")[[1]][3]
	tiss_use = strsplit(df_name, "_")[[1]][2]
	group_use = paste(sp_use, tiss_use, sep = "_")
	
	cons     = subset(df, df$cons_gene_by_tiss_ord == "YES")
	backgr   = subset(df, df$cons_gene_by_tiss_ord == "NO")
	
	cons_backgr_dnds  = t.test(cons$dnds, backgr$dnds,   paired = FALSE, alternative = c("two.sided"))
	cons_backgr_pnps  = t.test(cons$pnps, backgr$pnps,   paired = FALSE, alternative = c("two.sided"))

	
	dnds_t <- as.vector(c(cons_backgr_dnds$statistic))
	dnds_p <- c(cons_backgr_dnds$p.value)	
	pnps_t <- as.vector(c(cons_backgr_pnps$statistic))
	pnps_p <- c(cons_backgr_pnps$p.value)	
	
	group = rep(group_use, 1)
	comp = c("cons_backgr")
	
	df_out <- cbind(comp, group, dnds_t,dnds_p,pnps_t,pnps_p)
	return(df_out)
}

ttest_2_allTiss_overall <- as.data.frame(rbind(
ty_tests_2(dat1_WB),
ty_tests_2(dat1_RT),
ty_tests_2(dat1_LG)))

str(ttest_2_allTiss_overall)

ttest_2_allTiss_overall$dnds_t <- as.numeric(as.character(ttest_2_allTiss_overall$dnds_t))
ttest_2_allTiss_overall$pnps_t <- as.numeric(as.character(ttest_2_allTiss_overall$pnps_t))

WB_dnds_obs_tt_TS <- ttest_2_allTiss_overall$dnds_t[1]
RT_dnds_obs_tt_TS <- ttest_2_allTiss_overall$dnds_t[2]
LG_dnds_obs_tt_TS <- ttest_2_allTiss_overall$dnds_t[3]

WB_pnps_obs_tt_TS <- ttest_2_allTiss_overall$pnps_t[1]
RT_pnps_obs_tt_TS <- ttest_2_allTiss_overall$pnps_t[2]
LG_pnps_obs_tt_TS <- ttest_2_allTiss_overall$pnps_t[3]



get_adjusted_t_pval_dnds = function(rand_df,calc_TS){
	df_r = rand_df
	subset_larger = ""
	
	## two-tailed
	if(calc_TS > 0){
		df_r$dnds_TS_class = ifelse(df_r$dnds_t_teststat < calc_TS, "smaller", "larger")
		subset_larger = subset(df_r, df_r$dnds_TS_class == "larger")		
	}
	else{
		df_r$dnds_TS_class = ifelse(df_r$dnds_t_teststat < calc_TS, "larger", "smaller")
		subset_larger = subset(df_r, df_r$dnds_TS_class == "larger")		
	}

	## 
	
	tail_rank = length(subset_larger[,1])
	adj_pval = 10000
	if(tail_rank == 0){
		adj_pval = 0
	}
	else{
		adj_pval = tail_rank /10000 * 2
		
	}
	
	print(tail_rank )
	print(adj_pval)
	return(adj_pval)
}


get_adjusted_t_pval_pnps = function(rand_df,calc_TS){
	df_r = rand_df
	subset_larger = ""
	
	## two-tailed
	if(calc_TS > 0){
		df_r$pnps_TS_class = ifelse(df_r$pnps_t_teststat < calc_TS, "smaller", "larger")
		subset_larger = subset(df_r, df_r$pnps_TS_class == "larger")		
	}
	else{
		df_r$pnps_TS_class = ifelse(df_r$pnps_t_teststat < calc_TS, "larger", "smaller")
		subset_larger = subset(df_r, df_r$pnps_TS_class == "larger")		
	}

	## 
	
	tail_rank = length(subset_larger[,1])
	adj_pval = 10000
	if(tail_rank == 0){
		adj_pval = 0
	}
	else{
		adj_pval = tail_rank /10000 * 2
		
	}
	
	print(tail_rank )
	print(adj_pval)
	return(adj_pval)
}



dnds_perm_t_adj_pval <- c(
get_adjusted_t_pval_dnds(WB_t_rand_df, WB_dnds_obs_tt_TS),
get_adjusted_t_pval_dnds(RT_t_rand_df, RT_dnds_obs_tt_TS),
get_adjusted_t_pval_dnds(LG_t_rand_df, LG_dnds_obs_tt_TS))

pnps_perm_t_adj_pval <- c(
get_adjusted_t_pval_pnps(WB_t_rand_df, WB_pnps_obs_tt_TS),
get_adjusted_t_pval_pnps(RT_t_rand_df, RT_pnps_obs_tt_TS),
get_adjusted_t_pval_pnps(LG_t_rand_df, LG_pnps_obs_tt_TS))

ttest_allTiss_overall_adj_pval <- cbind(ttest_2_allTiss_overall,dnds_perm_t_adj_pval ,pnps_perm_t_adj_pval )
write.csv(ttest_allTiss_overall_adj_pval, "10sp_dndspnps_output/ttest_allTiss_overall_adj_pval.csv")


########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), "10sp_dndspnps_output/10sp_pNpSdNdS.R_sessionInfo.txt")






