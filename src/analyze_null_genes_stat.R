library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)



########################################### LOAD data ###################################
############ ** ASSOCIATED genes ** ##########
load("RData/analyze_null_genes_Broad_associated_priority.RData") #time_elapsed, list.par_analysis

############ ** PRIORITIZED genes ** ##########
#load("????")

############################# EXTRACT BROAD DATA #################################
############### Extracting from list
list.null.summary <- lapply(list.par_analysis, "[[", "df.null.summary")
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal <- lapply(list.null.fits, "[[", "fit.natal")
list.null.fit.prioritized.higher <- lapply(list.null.fits, "[[", "fit.prioritized.higher")
###### Combining
df.null.summary <- ldply(list.null.summary) # COMBINING list of data frames


############################# STATISTCAL TESTs #################################

################### ASSOCIATED genes #############
### *** REMEMBER to update the list.null.fit.natal before running the ttest

### Empirical p-value: NATAL TEST
tmp.obs.t_statistic <- -3.2976
tmp.null_distribution.t_statistic <- sapply(list.null.fit.natal, function(x) {x$statistic}) # same as sapply(list.null.natal_fits, "[[", c("statistic"))
sum(tmp.null_distribution.t_statistic > tmp.obs.t_statistic)/length(tmp.null_distribution.t_statistic) # calc empirical p-val, one sided, alternative="greater". WE ARE LOOKING AT RIGHT SIDED TAIL.
qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + coord_cartesian(ylim=c(0, 125))


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


