options(echo=TRUE)


library(foreach)
library(doMC)

library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext

rm(list=ls())

wd <- "/cvar/jhlab/timshel/scz"
setwd(wd)

###################################### Switch ######################################
################## FILES - ASSOCIATED VS PRIORITIZED ##################
path.local <- "/cvar/jhlab/timshel/scz"
file.null_genes <- paste(path.local, "/data/schizophrenia_expression325permutations0to999.genes.combined.csv", sep="") # Associated genes!
#file.null_genes <- paste(path.local, "/data/schizophrenia_expression325permutations0to999.genes.prioritized.top_63.combined.csv", sep="") # Prioritized genes!

################## Data ##################
### RNAseq
file.RData <- "RData/data_rnaseq_expression_processed.RData"
data_Rdata_result <- "null_RData_broad_rnaseq_associated_paired-ttest_priority.RData" # [associated or prioritized][priority or ish]
#data_Rdata_result <- "null_RData_broad_rnaseq_prioritized_paired-ttest_priority.RData" # prioritized

### Microarray
#file.RData <- "RData/data_marray_expression.RData"
#data_Rdata_result <- "null_RData_broad_marray_associated_paired-ttest_priority.RData" # [associated or prioritized][priority or ish]
#data_Rdata_result <- "null_RData_broad_marray_prioritized_paired-ttest_priority.RData" # prioritized


################## Print ##################
print(paste("data_Rdata_result:", data_Rdata_result))
print(paste("file.RData:", file.RData))


###################################### LOAD DATA ######################################
load(file.RData)

########################### NULL FILES ####################
######### Reading null file
df.null_genes <- read.csv(file.null_genes,h=T)
str(df.null_genes)
n_permutations <- length(unique(df.null_genes$permutation))


################ CPU Parameters #########################
print(paste("number of CPUs available:", detectCores()))
registerDoMC(4)
print(paste("number of CPUs used:", getDoParWorkers()))


####################### LOOP parameters #########################
foreach_min <- min(df.null_genes$permutation)
foreach_max <- max(df.null_genes$permutation)
#foreach_min <- 0
#foreach_max <- 6
print(paste("foreach_min:", foreach_min))
print(paste("foreach_max:", foreach_max))

### CONSIDER USING THE FOLLOWING CONSTRUCT INSTEAD [without the idx] (from analyze_null_genes_prioritized_tmp.R)
# idx = 0
# for (i in unique(df.null_genes$permutation)) {
# #for (i in 0:2) {
#   time.start <- proc.time()
#   idx <- idx+1
#   ....
#   df.null.current <- subset(df.null_genes, permutation==i)
#   df.null.mapping[idx,"n_mapped_genes"] <- ....

####################### Running PARALLEL loop #########################
time.start <- proc.time()
list.par_analysis <- foreach (i=foreach_min:foreach_max, .packages=c("plyr")) %dopar% {
  time.loop.start <- proc.time()
  i = i+1
  print(paste("i is:", i))

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
  
  str(df.expr.subset)
  
  ################## PAIRED t-test ##################
  df.natal.gene.mean <- ddply(df.expr.subset, .(ensembl_gene_id, natal), summarise,
                              mean=mean(value, na.rm=T))
  fit.natal.paired <- t.test(mean~natal, data=df.natal.gene.mean, alternative="greater", paired=TRUE)


  ######### T-test: prenatal vs postnatal ######
  ### *** DUMMY INSERT - UNCOMMENT ME LATER ** ###
  #fit.natal <- t.test(value~natal, data=df.expr.subset, alternative="greater")
  fit.natal <- "dummy_fit"

  ######### T-test: higher expressed ######
  ### *** DUMMY INSERT - UNCOMMENT ME LATER ** ###
  #df.expression_matrix.clean.melt$gene_type <- as.factor(ifelse(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.null.current$ensembl_gene_id, "prioritized", "other"))
  #fit.prioritized.higher <- t.test(value~gene_type, data=df.expression_matrix.clean.melt, alternative="less") # other < prioritized ("o" comes before "p")
  fit.prioritized.higher <- "dummy_fit"


  ######### T-test: higher expressed for each stage ######
  ### *** DUMMY INSERT - UNCOMMENT ME LATER ** ###
  # linmod <- function(df) {
  #   t.test(value~gene_type, data=df, alternative="less") # other < prioritized ("o" comes before "p")
  # }
  # list.fit.stage <- dlply(df.expression_matrix.clean.melt, .(stage), linmod) # this is a LIST of fits (one element per stage)
  list.fit.stage <- "dummy_fit"


  ########### Calculating summaries - **MEAN** ##########
  df.null.mean.summary <- ddply(df.expr.subset, c("stage", "structure_acronym"), summarise,
                      mean = mean(value, na.rm=TRUE),
                      sd   = sd(value, na.rm=TRUE))
  ### Adding factor
  df.null.mean.summary$permutation <- i

  ########### Calculating summaries - **MEDIAN** ##########
  df.null.median.summary <- ddply(df.expr.subset, c("stage", "structure_acronym"), summarise,
                      mean = median(value, na.rm=TRUE),
                      sd   = sd(value, na.rm=TRUE))
  ### Adding factor
  df.null.median.summary$permutation <- i
  
  ########################## Returning results ######################
  time.loop <- proc.time() - time.loop.start
  list.null.fits = list(fit.natal=fit.natal, fit.natal.paired=fit.natal.paired, fit.prioritized.higher=fit.prioritized.higher, list.fit.stage=list.fit.stage)
  list(df.null.mapping=df.null.mapping, df.null.mean.summary=df.null.mean.summary, df.null.median.summary=df.null.median.summary, list.null.fits=list.null.fits, time.loop=time.loop)
}
###### TIME for foreach loop
time_elapsed <- proc.time() - time.start
print("EXECUTION time for foreach loop:")
print(time_elapsed)
##### TIME mean for each loop
mean_loop_time <- mean(sapply(list.par_analysis, "[[", "time.loop"))
print(paste("MEAN loop time:", mean_loop_time))


#################################### Save to .Rdata ######################################
#save(time_elapsed, df.null.mapping, df.null.combined, list.null, list.null.natal_fits, file = data_Rdata_result)
print("Saving workspace #1")
save(time_elapsed, list.par_analysis, file=data_Rdata_result)

######################################## FINISH ###########################################
print("THE SCRIPT HAS COMPLETED")



