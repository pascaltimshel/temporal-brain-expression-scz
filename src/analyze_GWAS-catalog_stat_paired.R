library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)



########################################### LOAD data ###################################

############################# EXTRACT BROAD DATA #################################
source("function_GWAScatalog_get_summary.R")
################################ ** PRIORITIZED genes ** #########################
load("RData/GWAScatalog_prioritized_RData_broad_rnaseq_priority.RData")
names(list.par_analysis)
############### Extracting from list
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal.paired <- llply(list.null.fits, "[[", "fit.natal.paired") # *NEW PAIRED* 


############################# SUBSETTING on criteria #################################
####### csv file to filter on the number of genes - ONLY THE PRIORITIZED GENES SHOULD BE USED FOR SUBSETTING.
file.gwas_catalog_stat <- "../data/GWAS_catalog_prioritized_genes.csv"
#file.gwas_catalog_stat <- "../data/GWAS_catalog_associated_genes.csv"

####### read csv file
df.gwas_catalog_stat <- read.csv(file.gwas_catalog_stat)
str(df.gwas_catalog_stat)
df.gwas_catalog_stat.summary <- ddply(df.gwas_catalog_stat, .(permutation), summarise,
                                      n_genes = length(n_genes))
criteria <- 10
df.gwas_catalog_stat.summary.crit <- subset(df.gwas_catalog_stat.summary, n_genes>=criteria)
phenotypes_passing_criteria <- as.character(df.gwas_catalog_stat.summary.crit$permutation)



############################# STATISTCAL TESTs #################################
################### GWAS genes #############
### *** REMEMBER to update the list.null.fit.natal before running the ttest

list.null.fit.natal.paired.passing_criteria <- list.null.fit.natal.paired[names(list.null.fit.natal.paired) %in% phenotypes_passing_criteria]

### Empirical p-value: NATAL TEST
obs.t_statistic <- 2.6439 # Observed value for prioritized SCZ genes
df.null.t_statistic <- ldply(list.null.fit.natal.paired.passing_criteria, function(x) {t=x$statistic}, .id="id") # same as sapply(list.null.natal_fits, "[[", c("statistic"))
df.null.t_statistic
empirical.pvalue <- sum(df.null.t_statistic$t > obs.t_statistic)/length(df.null.t_statistic$t) # calc empirical p-val, one sided, alternative="greater"
empirical.pvalue
p <- ggplot(data=df.null.t_statistic, aes(x=t)) + geom_histogram(binwidth=1)
p <- p + geom_vline(xintercept=c(obs.t_statistic), linetype="dotted", size=1)
p <- p + labs(title=paste("GWAS Catalog. Empirical distribution. Criteria=", criteria, ": n_traits passing filter = ",nrow(df.null.t_statistic), "| pval=", round(empirical.pvalue,3), sep=""), y="Number of tests", x="t-statistic")
# + coord_cartesian(ylim=c(0, 7))
### text annotating
df.null.t_statistic.extreme <- subset(df.null.t_statistic, t>obs.t_statistic)
df.null.t_statistic.extreme$position <- df.null.t_statistic.extreme$t + runif(2,-0.5,0.5)
p + geom_text(data=df.null.t_statistic.extreme, aes(x=position, y=2.2, label=id), angle=60, hjust=0)



