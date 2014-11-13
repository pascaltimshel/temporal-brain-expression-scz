############### SYNOPSIS ###################
# This script is *DEPRECATED*. I am unsure if the script works. USE THE BROAD VERSION.
# The only reason to use this script is to test the PARELLEL for-loop locally (OSX).
############################################

library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext

### Parallel
library(foreach)
library(doMC)

rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

####################################### LOAD DATA ###########################################
###### LOADING EXPRESSION DATA
#source("function_read_marray.R", echo=TRUE)
# OUTPUT: df.expression_matrix.clean.melt + much more
# Loads data in "fresh" from text files
load(file="RData/data_marray_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
# gives the same is source("function_read_marray.R"), but loads data from RData


######################################## NULL FILES #######################################
######### Reading null file
file.null_genes <- "../Data/schizophrenia_expression325permutations0to999.genes.combined.csv"
#file.null_genes <- **** INSERT PATH TO PRIORITIZED genes *** #
df.null_genes <- read.csv(file.null_genes,h=T)

str(df.null_genes)
n_permutations <- length(unique(df.null_genes$permutation))


################ Running PARALLEL loop #########################
n_cpu <- detectCores()
print(n_cpu)
#registerDoMC(detectCores())
registerDoMC(n_cpu)
getDoParWorkers()

#list.par_analysis <- foreach (i=min(df.null_genes$permutation):max(df.null_genes$permutation)) %dopar% {
time_start <- proc.time()
list.par_analysis <- foreach (i=1:12) %dopar% {
  i = i+1
  ####### Subsetting null genes to current permutation #####
  df.null.current <- subset(df.null_genes, permutation==(i-1)) # OBS: i-1, because permutations are zero-based
  ##### Check mapping - compare to "df.expression_matrix.clean" NOT MOLTEN, will be too slow ####
  n_mapped_genes <- sum(df.expression_matrix.clean$ensembl_gene_id %in% df.null.current$ensembl_gene_id) # mapped genes
  n_unmapped_genes <- sum(!df.null.current$ensembl_gene_id %in% df.expression_matrix.clean$ensembl_gene_id) # non-mapped genes
  df.null.mapping <- data.frame(n_mapped_genes=n_mapped_genes, n_unmapped_genes=n_unmapped_genes)
  
  ########### Extracting expression data for current null genes ##########
  df.expr.subset <- subset(df.expression_matrix.clean.melt, subset=ensembl_gene_id %in% df.null.current$ensembl_gene_id)
  ##### adding stage_natal (prenatal vs post-natal) #####
  df.expr.subset$natal <- NA
  df.expr.subset[df.expr.subset$stage %in% c("s1","s2a","s2b","s3a","s3b","s4","s5"), "natal"] <- "prenatal"
  df.expr.subset[df.expr.subset$stage %in% c("s6","s7","s8","s9","s10","s11"), "natal"] <- "postnatal"
  df.expr.subset$natal <- as.factor(df.expr.subset$natal)
  
  ######### T-test: prenatal vs postnatal ######
  fit1 <- t.test(value~natal, data=df.expr.subset, alternative="greater")
  
  ########### Calculating summaries ##########
  df.summary <- ddply(df.expr.subset, c("stage", "structure_acronym"), summarise,
                      mean = mean(value, na.rm=TRUE),
                      sd   = sd(value, na.rm=TRUE))
  ### Adding factor
  df.summary$permutation <- i
  
  ########################## Returning results ######################
  list(df.null.mapping=df.null.mapping, list.null=df.summary, list.null.natal_fits=fit1)
}
time_elapsed <- proc.time() - time_start
print(time_elapsed)

########################### Check of result lists ###################################
names(list.par_analysis)

### Extracting from list
list.null <- lapply(list.par_analysis, "[[", "list.null")
list.null.natal_fits <- lapply(list.par_analysis, "[[", "list.null.natal_fits")
#df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])})
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping")
### Combining
df.null.combined <- ldply(list.null) # COMBINING list of data frames

### Save to .Rdata
save(time_elapsed, df.null.mapping, df.null.combined, list.null, list.null.natal_fits, file = "RData/analyze_null_genes_OSX_parallel.RData")
#### LOAD data
#load("analyze_null_genes.RData")
