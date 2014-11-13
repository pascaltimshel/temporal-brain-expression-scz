############### SYNAPSIS ###################
# Name: graphics_diagnostics_rnaseq.R
# This script will produce DENSITY PLOTs
# The script also investigates the WINSORIZATION and LOG2 transformation of the "raw" RNAseq data
############################################

library(plyr)
library(ggplot2)
library(reshape2)


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

################################# LOAD Data #######################################

############ RNAseq
load(file="RData/data_rnaseq_expression_unprocessed.RData") # df.expression_matrix.clean.unprocessed, df.expression_matrix.clean.melt.unprocessed
load(file="RData/data_rnaseq_expression_processed.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt

################################# DIAGNOSTICS - RNAseq #######################################
# ALL SAMPLES - OBS: ***this takes a long time for RNAseq data***
p <- ggplot(df.expression_matrix.clean.melt, aes(x=value, fill=variable)) + geom_density(alpha=.3)
p


### SINGLE SAMPLE density plot | 13058_Ocx_s2a
df.expression_matrix.clean.melt.subset <- subset(df.expression_matrix.clean.melt, variable=="13058_Ocx_s2a")
df.expression_matrix.clean.melt.subset <- df.expression_matrix.clean.melt.subset[order(-df.expression_matrix.clean.melt.subset$value),]

### Log2 density plot
p <- ggplot(df.expression_matrix.clean.melt.subset, aes(x=log2(value+1), fill=variable)) + geom_density(alpha=.3)
p

### "Winsorizing" samples higher than 50
sum(df.expression_matrix.clean.melt.subset$value > 50) # --> 1320
t <- quantile(df.expression_matrix.clean.melt.subset$value, probs=seq(0,1,0.01))
tmp <- df.expression_matrix.clean.melt.subset[df.expression_matrix.clean.melt.subset$value < 50,]
p <- ggplot(tmp, aes(x=log2(value+1), fill=variable)) + geom_density(alpha=.3)
p

