############### SYNOPSIS ###################
# This is a TEMPORARY script similar to "analyze_null_genes_mapping.R"
# Notice the loop structure can loop over "unfinished" permutations. That is, permutations where not all 1000 permutations are complete.
# 1) Plot the distribution of number of genes and mapped genes for NULL GWAS of prioritized genes
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)


######################################## NULL FILES #######################################
######### Reading null file
file.null_genes <- "/Users/pascaltimshel/p_scz/brainspan/src/schizophrenia_expression325permutations0to999.genes.prioritized.combined.csv"
df.null_genes <- read.csv(file.null_genes,h=T)
str(df.null_genes)
n_permutations <- length(unique(df.null_genes$permutation))
n_permutations

######## Exploratory analysis - unmapped genes
df.null.stats <- ddply(df.null_genes, c("permutation"), summarise,
      n = length(ensembl_gene_id),
      n_p005 = sum(pval_nominal<0.05))
### histogram plot
ggplot(df.null.stats, aes(x=n_p005)) + geom_histogram(aes(y = ..density..)) + geom_density() + labs(title="number of associated genes")
### range of associated loci
range(df.null.stats$n_p005)
mean(df.null.stats$n_p005) 


######################################## NULL FILES - top54 prioritized genes #######################################
######### Reading null file
file.null_genes <- "/Users/pascaltimshel/p_scz/brainspan/src/schizophrenia_expression325permutations0to999.genes.prioritized.top54.combined.csv"
df.null_genes <- read.csv(file.null_genes,h=T)
str(df.null_genes)
n_permutations <- length(unique(df.null_genes$permutation))
n_permutations

######## Exploratory analysis - unmapped genes
df.null.stats <- ddply(df.null_genes, c("permutation"), summarise,
                       n = length(ensembl_gene_id))
## check that all permutations have 54 genes


################################### FOR LOOP #############################
### Initializing data frame
df.null.mapping <- data.frame(n_mapped_genes=NA, n_unmapped_genes=NA)
### Running loop
idx = 0
for (i in unique(df.null_genes$permutation)) {
#for (i in 0:2) {
  time.start <- proc.time()
  idx <- idx+1
  #i <- i+1
  print(i)
  ####### Subsetting null genes to current permutation #####
  df.null.current <- subset(df.null_genes, permutation==i)
  ##### Check mapping - compare to "df.expression_matrix.clean" NOT MOLTEN, will be too slow ####
  df.null.mapping[idx,"n_mapped_genes"] <- sum(df.expression_matrix.clean$ensembl_gene_id %in% df.null.current$ensembl_gene_id) # mapped genes
  df.null.mapping[idx,"n_unmapped_genes"] <- sum(!df.null.current$ensembl_gene_id %in% df.expression_matrix.clean$ensembl_gene_id) # non-mapped genes
  #cat(as.character(df.null.current[!df.null.current$ensembl_gene_id %in% df.expression_matrix.clean$ensembl_gene_id, "ensembl_gene_id"]), sep="\n") # print non-mapped genes
}

########################################### ??? ###################################
######## Exploratory analysis - unmapped genes
#### histogram for mapped genes
ggplot(df.null.mapping, aes(x=n_mapped_genes)) + geom_histogram(aes(y = ..density..)) + geom_density() + labs(title="number of mapped prioritized genes (n_perm=28)")
###### SUMMARY stats for mapped genes
mean(df.null.mapping$n_mapped_genes) 
range(df.null.mapping$n_mapped_genes) 
sum(df.null.mapping$n_mapped_genes<47) 
