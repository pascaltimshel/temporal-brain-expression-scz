library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext



############################# EXTRACT BROAD DATA #################################

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
############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
############## Combining summary data frames
df.null.median.summary.prio <- ldply(list.null.median.summary) # COMBINING list of data frames
df.null.median.summary.sem.prio <- ddply(df.null.median.summary.prio, .(stage, permutation), summarise,
                                         mean1 = mean(mean, na.rm=TRUE),
                                         sd1   = sd(mean, na.rm=TRUE))



################################ ** ASSOCIATED genes ** #########################
load("RData/GWAScatalog_associated_RData_broad_rnaseq_priority.RData")
names(list.par_analysis)
############### Extracting from list
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal.paired <- llply(list.null.fits, "[[", "fit.natal.paired") # *NEW PAIRED* 
############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
############## Combining summary data frames
df.null.median.summary.assoc <- ldply(list.null.median.summary) # COMBINING list of data frames
df.null.median.summary.sem.assoc <- ddply(df.null.median.summary.assoc, .(stage, permutation), summarise,
                                          mean1 = mean(mean, na.rm=TRUE),
                                          sd1   = sd(mean, na.rm=TRUE))

var_keep <- c("wd", 
              "df.null.median.summary.prio", "df.null.median.summary.sem.prio",
              "df.null.median.summary.assoc", "df.null.median.summary.sem.assoc")
#list_rm <- setdiff(ls(), var_keep)
rm(list=setdiff(ls(), var_keep))

