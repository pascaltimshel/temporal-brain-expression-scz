############### SYNOPSIS ###################
# This script will make various plots for the NULL GWAS
# See the script "analyze_null_genes_mapping.R" for plotting distributions
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

########################################### LOAD data ###################################
############ ** ASSOCIATED genes ** ##########
load("RData/analyze_null_genes_Broad_associated_priority.RData") # list.par_analysis + more


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

######################################## NULL FILES #######################################

################## SUMMARIZE df.null.summary #######################

df.null.summary.stage <- ddply(df.null.summary, c("stage", "permutation"), summarise,
                               mean1 = mean(mean, na.rm=TRUE),
                               sd   = sd(mean, na.rm=TRUE))

df.null.summary.stage.mean <- ddply(df.null.summary, c("stage"), summarise,
                                    mean1 = mean(mean, na.rm=TRUE),
                                    sd   = sd(mean, na.rm=TRUE))



############################### PLOT - 1 ################################
p <- ggplot(df.null.summary.stage, aes(x=stage, y=mean1, group=permutation)) + geom_line(aes(colour = permutation))
#p <- p + geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd), width=0.2)
#p <- p + guides(colour=guide_legend(nrow = 10))
p <- p + labs(title = "permutation time series")
p
### Adding mean (base line)
p <- p + geom_line(data=df.null.summary.stage.mean, aes(x=stage, y=mean1, group=1), linetype="dashed", size=2.5)
p <- p + geom_errorbar(data=df.null.summary.stage.mean, aes(x=stage, group=1, ymax = mean1 + sd, ymin=mean1 - sd), width=0.2)
p

###### ADDING PRIORITIZED GENES (*******must run some code in df.summary.mean*****)
#p <- p + geom_line(data=df.summary.mean, aes(x=stage, y=mean, group = gene_type, linetype=gene_type), size=2.5)
#p

