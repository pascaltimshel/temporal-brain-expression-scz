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

################## Switch ##################
data_Rdata_result <- "analyze_null_genes_Broad_associated_priority.RData"
data_pre_initialization <- "pre_initialization_assocation_priority.RData"

print(paste("data_Rdata_result:", data_Rdata_result))
print(paste("data_pre_initialization:", data_pre_initialization))


run_initialization <- TRUE # IMPORTANT!
print(paste("Running initialization:", run_initialization))
############################################

if (run_initialization) {
  ################## PATHS ##################
  path.local <- "/cvar/jhlab/timshel/scz"
  file.expression_matrix <- paste(path.local, "/data/expression_matrix.csv", sep="")
  file.columns <- paste(path.local, "/data/columns_metadata.csv", sep="")
  file.rows <-  paste(path.local, "/data/rows_metadata.csv", sep="")
  ### IN src fold
  file.gene_length <- paste(path.local, "/data/df.clean.gene_length.insync.csv", sep="")
  file.null_genes <- paste(path.local, "/data/schizophrenia_expression325permutations0to999.genes.combined.csv", sep="")
  
  
  
  
  ################### Defining stages ####################
  stages = list()
  stages[["s1"]] = c() # 1 4-7 pcw Embryonic
  stages[["s2a"]] = c("8 pcw","9 pcw") # 2A 8-9 pcw Early prenatal
  stages[["s2b"]] = c("12 pcw") # 2B 10-12 pcw Early prenatal
  stages[["s3a"]] = c("13 pcw") # 3A 13-15 pcw Early mid-prenatal
  stages[["s3b"]] = c("16 pcw","17 pcw") # 3B 16-18 pcw Early mid-prenatal
  stages[["s4"]] = c("19 pcw","21 pcw","24 pcw") # 4 19-24 pcw Late mid-prenatal
  stages[["s5"]] = c("25 pcw","26 pcw","35 pcw","37 pcw") # 5 25-38 pcw Late prenatal
  stages[["s6"]] = c("4 mos") # 6 Birth-5 months Early infancy
  stages[["s7"]] = c("10 mos") # 7 6-18 months Late infancy
  stages[["s8"]] = c("1 yrs","2 yrs","3 yrs","4 yrs") # 8 19 months-5 yrs Early childhood
  stages[["s9"]] = c("8 yrs","11 yrs") # 9 6-11 yrs Late childhood
  stages[["s10"]] = c("13 yrs","15 yrs","18 yrs","19 yrs") # 10 12-19 yrs Adolescence
  stages[["s11"]] = c("21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs") # 11 20-60+ yrs Adulthood
  order.stages <- c("s1", "s2a", "s2b", "s3a", "s3b", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11")
  
  
  ########### READ columns file ###########
  df.columns <- read.csv(file.columns,h=T,row.names=1)
  #### Add stage column
  df.columns$stage <- as.factor(sapply(df.columns$age, function(x) {names(stages)[sapply(stages, function(stage) x %in% stage)]}))
  #### Sort factor levels of "stage"
  df.columns$stage <- factor(df.columns$stage, levels(df.columns$stage)[match(order.stages, levels(df.columns$stage))])
  
  ########### READ row file ###########
  df.rows <- read.csv(file.rows,h=T,row.names=1)
  
  ########### READ gene_length file ###########
  df.gene_length <- read.csv(file.gene_length,h=T)
  sum(is.na(df.gene_length$gene_length))
  
  ########### READ AND MANIPULATE expression file ###########
  ### ** THIS TAKES SOME TIME ** ###
  print("Loading in expression matrix")
  df.expression_matrix <- read.csv(file.expression_matrix,h=F,row.names=1) # HEADER FALSE
  ### Removing duplicates
  df.expression_matrix.clean <- df.expression_matrix[!(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE)), ]
  df.rows.clean <- df.rows[!(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE)), ]
  ### Setting column names - must be done first!
  colnames(df.expression_matrix.clean) <- with(df.columns, paste(donor_id, structure_acronym, stage, sep="_"))
  ### *** Normalizing expression matrix *** ###
  #df.expression_matrix.clean <- as.data.frame(scale(df.expression_matrix.clean)) # COLUMN NORMALIZATION
  #df.expression_matrix.clean <- (df.expression_matrix.clean-rowMeans(df.expression_matrix.clean))/apply(df.expression_matrix.clean,1,sd) # ROW NORMALIZATION
  ### *** NEW COLUMNS *** ###
  ### Setting ensemblID
  df.expression_matrix.clean$ensembl_gene_id <- df.rows.clean$ensembl_gene_id
  ### Setting priorizied factor
  #df.expression_matrix.clean$gene_type <- as.factor(ifelse(df.expression_matrix.clean$ensembl_gene_id %in% df.gene_prioritization[,1], "prioritized", "other"))
  ### Setting gene_length
  df.expression_matrix.clean$gene_length <- df.gene_length$gene_length
  df.expression_matrix.clean[15000,c("ensembl_gene_id","gene_length")] #---> must give gene_length=13522
  str(df.expression_matrix.clean,list.len=Inf)
  
  
  ############################### MANIPULAION - ALL GENES - full ###########################
  ### Melting dataframe
  df.expression_matrix.clean.melt <- melt(df.expression_matrix.clean, id=c("ensembl_gene_id", "gene_length"))
  #df.expression_matrix.clean.melt <- melt(df.expression_matrix.clean, id=c("ensembl_gene_id", "gene_type"))
  head(df.expression_matrix.clean.melt)
  
  ### Creating new variables from string
  print("Start variable spitting")
  variable_split <- strsplit(as.character(df.expression_matrix.clean.melt$variable), "_")
  df.expression_matrix.clean.melt$donor_id <- as.factor(sapply(variable_split, "[[", 1))
  df.expression_matrix.clean.melt$structure_acronym <- as.factor(sapply(variable_split, "[[", 2))
  df.expression_matrix.clean.melt$stage <- as.factor(sapply(variable_split, "[[", 3))
  ## sorting stage levels
  df.expression_matrix.clean.melt$stage <- with(df.expression_matrix.clean.melt, factor(stage, levels(stage)[match(order.stages, levels(stage))]))
  ## adding stage_natal (prenatal vs post-natal)
  df.expression_matrix.clean.melt$natal <- NA
  df.expression_matrix.clean.melt[df.expression_matrix.clean.melt$stage %in% c("s1","s2a","s2b","s3a","s3b","s4","s5"), "natal"] <- "prenatal"
  df.expression_matrix.clean.melt[df.expression_matrix.clean.melt$stage %in% c("s6","s7","s8","s9","s10","s11"), "natal"] <- "postnatal"
  df.expression_matrix.clean.melt$natal <- as.factor(df.expression_matrix.clean.melt$natal)
  str(df.expression_matrix.clean.melt)
  
  print("Saving workspace for PRE-INITIALIZATION")
  #save(list = ls(all = TRUE), file = data_pre_initialization)
  save.image(file = data_pre_initialization)
} else {
  #print(dir())
  #print(getwd())
  print("Loading workspace for PRE-INITIALIZATION")
  load(file=data_pre_initialization)
  #load(file="pre_initialization_assocation.RData")
}
#load(file="pre_initialization_assocation.RData")
###################################### Saving workspace ######################################



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
#foreach_max <- 2
print(paste("foreach_min:", foreach_min))
print(paste("foreach_max:", foreach_max))

####################### Running PARALLEL loop #########################
time.start <- proc.time()
list.par_analysis <- foreach (i=foreach_min:foreach_max, .packages=c("plyr")) %dopar% {
  time.loop.start <- proc.time()
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
  fit.natal <- t.test(value~natal, data=df.expr.subset, alternative="greater")
  
  ######### T-test: higher expressed ######
  df.expression_matrix.clean.melt$gene_type <- as.factor(ifelse(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.null.current$ensembl_gene_id, "prioritized", "other"))
  fit.prioritized.higher <- t.test(value~gene_type, data=df.expression_matrix.clean.melt, alternative="less") # other < prioritized ("o" comes before "p")
  
  str(df.expression_matrix.clean.melt)
  str(df.expr.subset)
  ########### Calculating summaries ##########
  df.null.summary <- ddply(df.expr.subset, c("stage", "structure_acronym"), summarise,
                           mean = mean(value, na.rm=TRUE),
                           sd   = sd(value, na.rm=TRUE))
  ### Adding factor
  df.null.summary$permutation <- i
  
  ########################## Returning results ######################
  time.loop <- proc.time() - time.loop.start
  list.null.fits = list(fit.natal=fit.natal, fit.prioritized.higher=fit.prioritized.higher)
  list(df.null.mapping=df.null.mapping, df.null.summary=df.null.summary, list.null.fits=list.null.fits, time.loop=time.loop)
}
###### TIME for foreach loop
time_elapsed <- proc.time() - time.start
print("EXECUTION time for foreach loop:")
print(time_elapsed)
##### TIME mean for each loop
mean_loop_time <- mean(sapply(list.par_analysis, "[[", "time.loop"))
print(paste("MEAN loop time:", mean_loop_time))

###### Saving once.. *******CONSIDER DELETING THIS! *********
print("Saving workspace #1")
save(time_elapsed, list.par_analysis, file=data_Rdata_result)

########################### Check of result lists ###################################
names(list.par_analysis)

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



################################################### Save to .Rdata ######################################
#save(time_elapsed, df.null.mapping, df.null.combined, list.null, list.null.natal_fits, file = data_Rdata_result)
print("Saving workspace #2")
save(time_elapsed, list.par_analysis, file=data_Rdata_result)
#### LOAD data
#load("analyze_null_genes.RData")
