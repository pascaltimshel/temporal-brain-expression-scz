library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)



########################################### LOAD data ###################################
############ ** ASSOCIATED genes ** ##########
### marray
#load("RData/null_RData_broad_marray_associated_priority.RData") #time_elapsed, list.par_analysis
### rnaseq
load("RData/null_RData_broad_rnaseq_associated_priority.RData") #time_elapsed, list.par_analysis

############ ** PRIORITIZED genes ** ##########
#load("????")

############################# EXTRACT BROAD DATA #################################
############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
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
#tmp.obs.t_statistic <- -3.2976 # marray
tmp.obs.t_statistic <- -6.6052 # RNAseq
tmp.null_distribution.t_statistic <- sapply(list.null.fit.natal, function(x) {x$statistic}) # same as sapply(list.null.natal_fits, "[[", c("statistic"))
sum(tmp.null_distribution.t_statistic > tmp.obs.t_statistic)/length(tmp.null_distribution.t_statistic) # calc empirical p-val, one sided, alternative="greater". WE ARE LOOKING AT RIGHT SIDED TAIL.
qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + coord_cartesian(ylim=c(0, 125)) + labs(title="associated genes - one-sided; consider right tail")


### Empirical p-value: HIGHER EXPRESSION
tmp.obs.t_statistic <- -126.5334
tmp.null_distribution.t_statistic <- sapply(list.null.fit.prioritized.higher, function(x) {x$statistic}) # same as sapply(list.null.natal_fits, "[[", c("statistic"))
sum(tmp.null_distribution.t_statistic < tmp.obs.t_statistic)/length(tmp.null_distribution.t_statistic) # calc empirical p-val, one sided, alternative="less". WE ARE LOOKING AT LEFT SIDED TAIL.
qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + coord_cartesian(ylim=c(0, 125))



################### PRIORITIZED genes #############
### *** REMEMBER to update the list.null.fit.natal before running the ttest

### Empirical p-value: NATAL TEST
tmp.obs.t_statistic <- 18.8815
tmp.null_distribution.t_statistic <- sapply(list.null.fit.natal, function(x) {x$statistic}) # same as sapply(list.null.natal_fits, "[[", c("statistic"))
sum(tmp.null_distribution.t_statistic > tmp.obs.t_statistic)/length(tmp.null_distribution.t_statistic) # calc empirical p-val, one sided, alternative="greater"
qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + coord_cartesian(ylim=c(0, 125))


### Empirical p-value: HIGHER EXPRESSION
tmp.obs.t_statistic <- -142.2306
tmp.null_distribution.t_statistic <- sapply(list.null.fit.prioritized.higher, function(x) {x$statistic}) # same as sapply(list.null.natal_fits, "[[", c("statistic"))
sum(tmp.null_distribution.t_statistic < tmp.obs.t_statistic)/length(tmp.null_distribution.t_statistic) # calc empirical p-val, one sided, alternative="less". WE ARE LOOKING AT LEFT SIDED TAIL.
qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + coord_cartesian(ylim=c(0, 125))


