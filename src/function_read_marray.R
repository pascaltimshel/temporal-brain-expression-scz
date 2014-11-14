############### SYNOPSIS ###################
# Function to read MICROARRAY expression file
# This function similar to "function_read_rnaseq.R" expect for the INPUT and OUTPUT (.RData) files
# READ MICROARRAY: the microarray expression data does not need normalization. However, the ENSEMBL gene identifiers are CLEANED (duplicates removed)
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext

########### SOURCING ###########
source("function_def_stages.R", echo=TRUE)

########### FILES ##############
file.expression_matrix <- "../data/141031/microarray/expression_matrix.csv"
file.columns <- "../data/141031/microarray/columns_metadata.csv"
file.rows <- "../data/141031/microarray/rows_metadata.csv"
### IN src fold
file.gene_length <- "../data/annotation.gene_length.microarray_expression_df_clean.synced.csv"

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
df.expression_matrix <- read.csv(file.expression_matrix,h=F,row.names=1) # HEADER FALSE
### Removing duplicates
df.expression_matrix.clean <- df.expression_matrix[!(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE)), ]
df.rows.clean <- df.rows[!(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE)), ]
# Print the sum of duplicates
n_duplicates <- sum(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE))
print(paste("Number of ENSEMBL Gene IDs (duplicates) removed:", n_duplicates))
### Setting column names - must be done first!
colnames(df.expression_matrix.clean) <- with(df.columns, paste(donor_id, structure_acronym, stage, sep="_"))

### *** Normalizing expression matrix *** ###
#df.expression_matrix.clean <- as.data.frame(scale(df.expression_matrix.clean)) # COLUMN NORMALIZATION
#df.expression_matrix.clean <- (df.expression_matrix.clean-rowMeans(df.expression_matrix.clean))/apply(df.expression_matrix.clean,1,sd) # ROW NORMALIZATION

### *** NEW COLUMNS *** ###
### Setting ensemblID
df.expression_matrix.clean$ensembl_gene_id <- df.rows.clean$ensembl_gene_id

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

################################# SAVING DATA ###################################
#save(df.expression_matrix.clean, df.expression_matrix.clean.melt, file="RData/data_marray_expression.RData")
