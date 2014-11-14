############### SYNOPSIS ###################
# This script is *DEPRECATED* but works.
# The only reason to use this script is to test the for-loop locally (OSX).
# The script is useful for producing data for inspecting gene mapping distirbutions for NULL GWAS.
############################################

library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


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
file.null_genes <- "../data/schizophrenia_expression325permutations0to999.genes.combined.csv"
#file.null_genes <- **** INSERT PATH TO PRIORITIZED genes *** #
df.null_genes <- read.csv(file.null_genes,h=T)


################################### FOR LOOP #############################
### Initializing data frame
df.null.mapping <- data.frame(n_mapped_genes=NA, n_unmapped_genes=NA)
#df.null.mean <- data.frame(matrix(NA, nrow = n_permutations, ncol = 12))
list.null <- list()
list.null.natal_fits <- list()
### Running loop
for (i in min(df.null_genes$permutation):max(df.null_genes$permutation)) {
  #for (i in 0:2) {
  time.start <- proc.time()
  i = i+1
  print(i)
  ####### Subsetting null genes to current permutation #####
  df.null.current <- subset(df.null_genes, permutation==(i-1)) # OBS: i-1, because permutations are zero-based
  ##### Check mapping - compare to "df.expression_matrix.clean" NOT MOLTEN, will be too slow ####
  df.null.mapping[i,"n_mapped_genes"] <- sum(df.expression_matrix.clean$ensembl_gene_id %in% df.null.current$ensembl_gene_id) # mapped genes
  df.null.mapping[i,"n_unmapped_genes"] <- sum(!df.null.current$ensembl_gene_id %in% df.expression_matrix.clean$ensembl_gene_id) # non-mapped genes
  #cat(as.character(df.null.current[!df.null.current$ensembl_gene_id %in% df.expression_matrix.clean$ensembl_gene_id, "ensembl_gene_id"]), sep="\n") # print non-mapped genes
  
  ########### Extracting expression data for current null genes ##########
  #df.expr.subset <- subset(df.expression_matrix.clean, subset=ensembl_gene_id %in% df.null.current$ensembl_gene_id)
  df.expr.subset <- subset(df.expression_matrix.clean.melt, subset=ensembl_gene_id %in% df.null.current$ensembl_gene_id)
  #print(dim(df.expr.subset))
  ##### adding stage_natal (prenatal vs post-natal) #####
  df.expr.subset$natal <- NA
  df.expr.subset[df.expr.subset$stage %in% c("s1","s2a","s2b","s3a","s3b","s4","s5"), "natal"] <- "prenatal"
  df.expr.subset[df.expr.subset$stage %in% c("s6","s7","s8","s9","s10","s11"), "natal"] <- "postnatal"
  df.expr.subset$natal <- as.factor(df.expr.subset$natal)
  ######### T-test: prenatal vs postnatal ######
  ## ** UNCOMMENT THIS LATER ** 
  #fit1 <- t.test(value~natal, data=df.expr.subset, alternative="greater")
  fit1 <- "dummy - time saver" # 
  list.null.natal_fits[[i]] <- fit1
  
  ########### Calculating summaries ##########
  ## ** UNCOMMENT THIS LATER ** 
#   df.summary <- ddply(df.expr.subset, c("stage", "structure_acronym"), summarise,
#                       mean = mean(value, na.rm=TRUE),
#                       sd   = sd(value, na.rm=TRUE))
  df.summary <- data.frame(dummy1=rnorm(10), dummy2=rnorm(10))
  
  ### Adding factor
  df.summary$permutation <- i
  list.null[[i]] <- df.summary
  names(list.null)[i] <- paste("perm",i,sep="")
  time.elapsed <- proc.time() - time.start
  print(time.elapsed)
}

######## COMBINING list of data frames
df.null.combined <- ldply(list.null)


####################################### Save to .Rdata ##################################
#save(df.null.combined, list.null, list.null.natal_fits, file = "analyze_null_genes.RData")
save.image(file="analyze_null_genes_OSX.RData")
