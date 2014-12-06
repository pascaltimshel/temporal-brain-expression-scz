library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

source("multiplot.R")

########################################### LOAD data ###################################
############ ** ASSOCIATED genes ** ##########
### rnaseq
#load("RData/null_RData_broad_rnaseq_associated_paired-ttest_priority_ish.RData") #time_elapsed, list.par_analysis

############ ** PRIORITIZED genes ** ##########
### rnaseq
load("RData/null_RData_broad_rnaseq_prioritized_paired-ttest_priority_ish.RData") #time_elapsed, list.par_analysis


############################# EXTRACT BROAD DATA #################################
############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal.paired <- lapply(list.null.fits, "[[", "fit.natal.paired") # *NEW PAIRED* 
list.null.fit.natal <- lapply(list.null.fits, "[[", "fit.natal")
list.null.fit.prioritized.higher <- lapply(list.null.fits, "[[", "fit.prioritized.higher")
list.null.list.fit.stage <- lapply(list.null.fits, "[[", "list.fit.stage")
############## Extracting from list.null.list.fit.stage
df.fit.stage<-ldply(seq_along(list.null.list.fit.stage), function(i) {ldply(list.null.list.fit.stage[[i]], function(fit.stage) {data.frame(t_statistic=fit.stage$statistic, perm=i)})})

############## Combining summary data frames
df.null.mean.summary <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary <- ldply(list.null.median.summary) # COMBINING list of data frames


############################# STATISTCAL TESTs #################################

################### ASSOCIATED genes #############
### *** REMEMBER to update the list.null.fit.natal before running the ttest

### Empirical p-value: NATAL TEST
obs.t_statistic <- -1.3412 # SCZ associated genes, RNAseq
df.null.t_statistic <- ldply(list.null.fit.natal.paired, function(x) {t=x$statistic}, .id="id") # same as sapply(list.null.natal_fits, "[[", c("statistic"))
df.null.t_statistic
empirical.pvalue <- sum(df.null.t_statistic$t > obs.t_statistic)/length(df.null.t_statistic$t) # calc empirical p-val, one sided, alternative="greater"
empirical.pvalue
p <- ggplot(data=df.null.t_statistic, aes(x=t)) + geom_histogram() #binwidth=1
p <- p + geom_vline(xintercept=c(obs.t_statistic), linetype="dotted", size=1)
p <- p + labs(title=paste("Null GWAS - associated genes. pval=", round(empirical.pvalue,3), sep=""), y="Number of tests", x="t-statistic")
p
# + coord_cartesian(ylim=c(0, 7))


################### PRIORITIZED genes #############
### *** REMEMBER to update the list.null.fit.natal before running the ttest


### Empirical p-value: NATAL TEST
obs.t_statistic <- 2.6439 # Observed value for prioritized SCZ genes
df.null.t_statistic <- ldply(list.null.fit.natal.paired, function(x) {t=x$statistic}, .id="id") # same as sapply(list.null.natal_fits, "[[", c("statistic"))
df.null.t_statistic
empirical.pvalue <- sum(df.null.t_statistic$t > obs.t_statistic)/length(df.null.t_statistic$t) # calc empirical p-val, one sided, alternative="greater"
empirical.pvalue
p <- ggplot(data=df.null.t_statistic, aes(x=t)) + geom_histogram(binwidth=1)
p <- p + geom_vline(xintercept=c(obs.t_statistic), linetype="dotted", size=1)
p <- p + labs(title=paste("Null GWAS - prioritized genes. pval=", round(empirical.pvalue,3), sep=""), y="Number of tests", x="t-statistic")
# + coord_cartesian(ylim=c(0, 7))


#multiplot(p, p1, cols=2)
#empirical_distribution_null_GWAS_rnaseq_MULTIPLOT{1x2}{prioANDassoc}_genes_paired_test-10x6.pdf
#empirical_distribution_null_GWAS_rnaseq_associated_genes_paired_test-8x6.pdf
#empirical_distribution_rnaseq_prioritized_genes_prenatal-vs-postnatal-10x6


##################################### PLOTTING ######################################
########### LOAD EXPRESSION DATA + GENE LISTS DATA
load(file="RData/data_rnaseq_expression_processed_full_res.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt
str(df.expression_matrix.clean.melt)
source("function_load_gene_list_data_individual_gene_plots.R", echo=TRUE)
# --> do_stage_converter
# --> df.summary
# --> df.summary.sem
# --> df.all.sem

########## SUMMARIZE NULL GWAS DATA - prepare for plotting
# Stage average across structures FOR EACH PERMUTATION. 
df.null.median.summary.per_perm <- ddply(df.null.median.summary, .(permutation, stage), summarise,
                                         mean_stage_across_structure = mean(mean, na.rm=TRUE),
                                         sd_stage_across_structure   = sd(mean, na.rm=TRUE))
str(df.null.median.summary.per_perm)
# "Grand mean" for all permutations - *THIS IS CURRENTLY NOT USED*
df.null.median.summary.sem <- ddply(df.null.median.summary.per_perm, .(stage), summarise,
                                    mean1 = mean(mean_stage_across_structure, na.rm=TRUE),
                                    sd1   = sd(mean_stage_across_structure, na.rm=TRUE))

### Subset on more extreme t-statistics
df.null.t_statistic.extreme <- subset(df.null.t_statistic, t>obs.t_statistic)
df.null.t_statistic.extreme$permutation <- as.numeric(rownames(df.null.t_statistic.extreme))
str(df.null.t_statistic.extreme)

df.null.median.summary.per_perm.extreme <- subset(df.null.median.summary.per_perm, permutation %in% df.null.t_statistic.extreme$permutation)



p <- ggplot()
### Plot each EXTREME null GWAS
p <- p + geom_line(data=df.null.median.summary.per_perm.extreme, aes(x=stage, y=mean_stage_across_structure, group=permutation, color="Null GWAS with more\nextreme t-statistic"), linetype='solid', size=1)
p
### Adding mean Prioritized
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#d7191c', width=0.2)
p
### Adding mean ALL (df.all.sem)
p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
p

cmap_color <- c("Prioritized genes"="#d7191c",
                "All genes"="black",
                "Null GWAS with more\nextreme t-statistic"="#00BF7D")
p <- p + scale_color_manual(name="Gene list", values=cmap_color)
p
#p <- p + guides(color=guide_legend(title=""))
#rnaseq_plot_null_gwas_with_extreme_t-statistic-8x6

################################### HELPER FUNCTIONS ###############################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

cols = gg_color_hue(10)
cols
scales:::show_col(cols)

