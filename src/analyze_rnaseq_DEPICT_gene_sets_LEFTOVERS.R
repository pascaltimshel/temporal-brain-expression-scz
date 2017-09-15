############### SYNOPSIS ###################
# Extracted subpart of "analyze_rnaseq_DEPICT_gene_sets.R"
# *OWN FUNCTION WAS NEVER COMPLETED!*

#####################################################################################
######################### Correlate Genes - OWN FUNCTION ###########################
#####################################################################################
### **write code to test for ONE gene_set!!!***

###### Correlation per stage ######## 
### initialization
df.res <- data.frame()
#bind_rows(df.res, data.frame(y = 1:4))
### DO SOMETHING HERE!

df.brainspan <- df.summary.stage
for (i.stage in levels(df.brainspan$stage)) {
  print(i.stage)
  ### Subset
  df.brainspan.sub <- subset(df.brainspan, stage==i.stage)
  ### Match
  idx_order <- match(df.brainspan.sub$ensembl_gene_id, df.depict.gene_set$ensembl_gene_id)
  idx_order
  stopifnot(all(idx_order)) # MUST be TRUE: then we know that all genes have a match
  ### *ALIGNING* ENSEMBL gene IDs
  df.brainspan.sub <- df.brainspan.sub[idx_order,]
  stopifnot(all(df.brainspan.sub$ensembl_gene_id == df.depict.gene_set$ensembl_gene_id)) # MUST be TRUE: then we know all genes are ALIGNED
  
  ### Looping over GS
  cols.gs <- colnames(df.depict.gene_set)[!colnames(df.depict.gene_set) %in% "ensembl_gene_id"]
  for (i.gs in cols.gs) {
    
  }
}

str(df.brainspan.sub$ensembl_gene_id)
##########################################






#####################################################################################
################################### CODE ARCHIVE ####################################
#####################################################################################

# weighted mean: weighted.mean()

#####################################################################################
########################### Correlation: basic info ############################
#####################################################################################

?cor
?cor.test

x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
y <- c( 2.6,  3.1,  NA,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)
cor.test(x, y, na.action="na.omit") # na.action="na.omit" --> will only keep PAIRWISE COMPLETE OBSERVATIONS


#####################################################################################
########################### Join df.summary.*; melt afterwards ############################
#####################################################################################
### Join
system.time(df.summary <- full_join(df.summary.stage, df.summary.structure_acronym, by="ensembl_gene_id"))
# Note that the "mean" columns in the output are disambiguated with a suffix ("x", and "y")
# Runtime: ~30 sec?
df.summary
# Source: local data frame [6,235,944 x 5]
# Groups: ensembl_gene_id

#    ensembl_gene_id stage   mean.x structure_acronym   mean.y
# 1  ENSG00000000003   s2a 4.640608               A1C 2.619151
# 2  ENSG00000000003   s2a 4.640608               AMY 3.114684

### Melt --> *THIS IS NOT HELPFUL*. "variable" will contain factor levels of "mean.x" and "mean.y"
df.summary.melt <- melt(df.summary, id.vars=c("ensembl_gene_id", "stage", "structure_acronym"))
head(df.summary.melt)
# --> 12,471,888 obs. of  5 variables
# ensembl_gene_id stage structure_acronym variable    value
# 1 ENSG00000000003   s2a               A1C   mean.x 4.640608
# 2 ENSG00000000003   s2a               AMY   mean.x 4.640608


#####################################################################################
################################### Joining data ####################################
#####################################################################################

full_join(x, y) #--> should give two new columns with all
inner_join(x, y) # --> same as full_join(). This is because the two data frames contain exactly the same ENSEMBL gene IDs. If df.x had more (or less) gene IDs, the resulting joined df would be smaller.
# --> **TRY to do a inner_join with "df.expression_matrix.clean.melt" instead. Should give same result because intersection is used.
semi_join(x, y)

# If a match is not unique, a join will add all possible combinations (the Cartesian product) of the matching observations:
# ---> use semi_join() to avoid duplications; semi_join() only ever remove observations.

### Full join
system.time(df.summary.stage.full_join <- df.summary.stage %>% full_join(df.depict.gene_set.melt, by="ensembl_gene_id"))
# --> 234.917 s 
# Source: local data frame [34,297,692 x 5]
# Groups: ensembl_gene_id
# 
#    ensembl_gene_id stage     mean                      gene_set gene_set_z_score
# 1  ENSG00000000003   s2a 4.640608   Decreased.Vertical.Activity         1.985478
# 2  ENSG00000000003   s2a 4.640608                      Dendrite        -1.119173
# 3  ENSG00000000003   s2a 4.640608 Abnormal.Locomotor.Activation         1.444252

### Inner join
df.summary.stage.inner_join <- df.summary.stage %>% inner_join(df.depict.gene_set.melt, by="ensembl_gene_id")
# Source: local data frame [34,297,692 x 5]
# Groups: ensembl_gene_id
# 
#    ensembl_gene_id stage     mean                      gene_set gene_set_z_score
# 1  ENSG00000000003   s2a 4.640608   Decreased.Vertical.Activity         1.985478
# 2  ENSG00000000003   s2a 4.640608                      Dendrite        -1.119173

### Semi_join
system.time(df.summary.stage.semi_join <- df.summary.stage %>% semi_join(df.depict.gene_set.melt, by="ensembl_gene_id"))
# --> 77.117 s
# Source: local data frame [239,844 x 3]
# Groups: ensembl_gene_id
# 
# ensembl_gene_id stage     mean
# 1  ENSG00000000003   s2a 4.640608
# 2  ENSG00000000003   s2b 3.098992




#####################################################################################
############################ align_and_correlate_genes ##################################
align_and_correlate_genes <- function(df.brainspan, df.depict.gs) {
  ### Input ###
  # df.brainspan:   BrainSpan (subset) *MOLTEN and SUMMERIZED* data frame. [Long format]
  # df.depict.gs:   DEPICT GS matrix. Dimensions = [Genes X GS]
  ### OBS ###
  
  return(1)
}
x <- align_and_correlate_genes(df.brainspan=df.brainspan.sub, df.depict.gs=df.depict.gene_set)



#####################################################################################
########### Subset genes in "df.expression_matrix.clean" to those in DEPICT #########
head(df.expression_matrix.clean)

### Subsetting data frame
df.expression_matrix.clean.sub.depict <- subset(df.expression_matrix.clean, ensembl_gene_id %in% df.depict.gene_set[,1])
### Converting "ensembl_gene_id" to *CHARACTER* 
# UNSURE about the following statement --> this is important for correct reordering AND checking the all(.) statement
df.expression_matrix.clean.sub.depict$ensembl_gene_id <- as.character(df.expression_matrix.clean.sub.depict$ensembl_gene_id)

### Finding matches
idx_order <- match(df.expression_matrix.clean.sub.depict$ensembl_gene_id, df.depict.gene_set[,1])
idx_order
all(idx_order) # MUST be TRUE: then we know that all genes have a match
### *ALIGNING* ENSEMBL gene IDs
df.expression_matrix.clean.sub.depict <- df.expression_matrix.clean.sub.depict[idx_order,]
all(df.expression_matrix.clean.sub.depict$ensembl_gene_id == df.depict.gene_set[,1]) # MUST be TRUE: then we know all genes are ALIGNED


#####################################################################################
###### *ddply* - Median per stage FOR EACH GENE ######## 
library(foreach)
library(doMC)
print(paste("number of CPUs available:", detectCores()))
registerDoMC(4)
print(paste("number of CPUs used:", getDoParWorkers()))
system.time(df.summary.stage <- ddply(df.expression_matrix.clean.melt.sub.depict, .(ensembl_gene_id, stage), summarise,
                                      mean = median(value, na.rm=TRUE),
                                      sd   = sd(value, na.rm=TRUE), .parallel=TRUE))
# Time (non-parallel) --> 110.975 (elapsed)
# Time (8 x parallel) --> 272.940 (elapsed), 533 (user)

###### Notes for dply ######## 
# Nested ddply
### REMEMBER: the order of splitting DOES NOT matter. The data is split into all possible factor level combinations and the result is computed
#.(ensembl_gene_id, structure_acronym)
#.(ensembl_gene_id, stage, structure_acronym)
### REMEMBER: we cannot rely on ddply/dply to keep the correct ordering of the genes. That is, GOOD PRACTICE WOULD BE TO MATCH GENES "MANUALLY" INSIDE THE FUNCTION DDPLY RUNS when correlating genes.

#####################################################################################
#####################################################################################









