############### SYNOPSIS ###################
# This script will plot the distribution of number of genes and mapped genes for NULL GWAS
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
.......
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") 


######################################## NULL FILES #######################################
######### Reading null file
file.null_genes <- "../data/schizophrenia_expression325permutations0to999.genes.combined.csv"
#file.null_genes <- "../data/schizophrenia_expression325permutations0to999.genes.prioritized.top54.combined.csv" # PATH TO PRIORITIZED genes *** #

df.null_genes <- read.csv(file.null_genes,h=T)

######### Checking null file
str(df.null_genes)
n_permutations <- length(unique(df.null_genes$permutation))

######## Exploratory analysis - NUMBER OF GENES IN NULL GWAS
df.null.stats <- ddply(df.null_genes, c("permutation"), summarise,
                       n = length(ensembl_gene_id))
### histogram plot
ggplot(df.null.stats, aes(x=n)) + geom_histogram(aes(y = ..density..)) + geom_density() + labs(title="number of associated genes")
### range of associated loci
range(df.null.stats$n) # range: 348 492
mean(df.null.stats$n) # mean: 398.31
sum(df.null.stats$n<363) # number of permutations with less than 363 associations: 19

############################ *** After running for loop *** ##########################
######## Exploratory analysis - NUMBER OF MAPPED GENES FOR NULL GWAS
#### histogram for mapped genes
ggplot(df.null.mapping, aes(x=n_mapped_genes)) + geom_histogram(aes(y = ..density..)) + geom_density() + labs(title="number of mapped associated genes")
###### SUMMARY stats for mapped genes
mean(df.null.mapping$n_mapped_genes) # mean: 302.236
range(df.null.mapping$n_mapped_genes) # range: 251 383
sum(df.null.mapping$n_mapped_genes<290) # number of permutations with less than 290 associations: 242


