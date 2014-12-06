############### SYNAPSIS ###################
# This script is the EXTENDED VERSION of "graphics_main1_median.R"
# The script loads the same gene lists, but also load the "null" 
# The script will make additional statistical tests to plot
# Loads data for either RNAseq/Microarray and associated/prioritized null genes
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

############################# FUNCTIONS #################################

######### Adding x-tickmarks for stage
do_stage_converter <- function (p) {
  stage_converter <- c("s1"="Embryonic",
                       "s2a"="Early prenatal",
                       "s2b"="Early prenatal",
                       "s3a"="Early mid-prenatal",
                       "s3b"="Early mid-prenatal",
                       "s4"="Late mid-prenatal",
                       "s5"="Late prenatal",
                       "s6"="Early infancy",
                       "s7"="Late infancy",
                       "s8"="Early childhood",
                       "s9"="Late childhood",
                       "s10"="Adolescence",
                       "s11"="Adulthood")
  p <- p + scale_x_discrete(name="", labels = stage_converter) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=rel(1.15)))
  return(p)
}


############################# LOAD EXPRESSION DATA #################################
#load(file="RData/data_marray_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
load(file="RData/data_rnaseq_expression_processed_high_res.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt
str(df.expression_matrix.clean.melt)

######### *** NEW - 12/01/2014 ** ######## Defining new stages
# source high res stages
source("function_def_stages_high_res.R", echo=TRUE)

## sorting AGE levels
df.expression_matrix.clean.melt$age <- with(df.expression_matrix.clean.melt, factor(age, levels(age)[match(order.age, levels(age))]))
levels(df.expression_matrix.clean.melt$age)

########################################### LOAD NULL data ###################################
################################## ** ASSOCIATED genes ** ##########################
### marray
#load("RData/null_RData_broad_marray_associated_priority.RData") #time_elapsed, list.par_analysis
### rnaseq
load("RData/null_RData_broad_rnaseq_associated_priority.RData") #time_elapsed, list.par_analysis

############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
############## Combining summary data frames
df.null.mean.summary.assoc <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary.assoc <- ldply(list.null.median.summary) # COMBINING list of data frames
#*OBS*: consider doing one more ddply around a c("stage", "permutation")
df.null.median.summary.sem.assoc <- ddply(df.null.median.summary.assoc, c("stage"), summarise,
                                          mean1 = mean(mean, na.rm=TRUE),
                                          sd1   = sd(mean, na.rm=TRUE))

################################ ** PRIORITIZED genes ** #########################
### marray
#load("RData/null_RData_broad_marray_prioritized_priority.RData") #time_elapsed, list.par_analysis
### rnaseq
load("RData/null_RData_broad_rnaseq_prioritized_priority.RData") #time_elapsed, list.par_analysis

############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
############## Combining summary data frames
df.null.mean.summary.prio <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary.prio <- ldply(list.null.median.summary) # COMBINING list of data frames
df.null.median.summary.sem.prio <- ddply(df.null.median.summary.prio, c("stage"), summarise,
                                         mean1 = mean(mean, na.rm=TRUE),
                                         sd1   = sd(mean, na.rm=TRUE))

############################# READING GENE LISTs #################################
path.datafiles <- '/Users/pascaltimshel/p_scz/brainspan/gene_lists'

###### Read into a list of files - PATTERN VERSION - read ALL .txt files in directory:
#files <- list.files(path = path.datafiles, pattern = "*.txt", full.names = TRUE) #full path
#names(files) <- list.files(path = path.datafiles, pattern = "*.txt") # filename
#cat(names(files), sep="\n")

###### Read SPECIFIC FILES:
filenames2read <- c("gene_associated.txt", "gene_nearest.txt", "gene_prioritization.txt", "gene_psd_human.txt", "gene_psd_mouse.txt")
filenames2read <- c(filenames2read, "gilman_nn_2012_cluster1.ens", "gilman_nn_2012_cluster1a.ens", "gilman_nn_2012_cluster1b.ens", "gilman_nn_2012_cluster2.ens")
filenames2read <- c(filenames2read, "gulsuner_S3A_damaging_cases.ens")
files <- as.list(paste(path.datafiles, filenames2read, sep="/"))
names(files) <- filenames2read
files
list_of_data <- llply(files, read.csv)#row.names = 1 --> NO!, stringsAsFactors = FALSE

names(list_of_data)

extract_genes_from_molten_df <- function(df_gene_list) {
  print("done")
  df <- subset(df.expression_matrix.clean.melt, ensembl_gene_id %in% df_gene_list[,1])
}
df.gene_list <- ldply(list_of_data, extract_genes_from_molten_df, .id="gene_list")
## Converting .id=gene_list to factor
df.gene_list$gene_list <- as.factor(df.gene_list$gene_list) 
str(df.gene_list)
levels(df.gene_list$gene_list)

###################################### PROCESSING GILMAN ################################

gilman_lvl <- c("gilman_nn_2012_cluster1.ens", "gilman_nn_2012_cluster2.ens") # THIS should correspond to the filesnames

###### Mean per stage/structure/cluster
df.gilman.summary.cluster.sem <- ddply(ddply(subset(df.gene_list, gene_list %in% gilman_lvl), .(stage, structure_acronym, gene_list), summarise, mean=median(value, na.rm=TRUE)), .(stage, gene_list), summarise, mean1=mean(mean, na.rm=TRUE),  sd1=sd(mean, na.rm=TRUE))
## plyr magic for renaming column and factor level
df.gilman.summary.cluster.sem <- rename(df.gilman.summary.cluster.sem, c("gene_list"="cluster")) # column
df.gilman.summary.cluster.sem$cluster <- revalue(df.gilman.summary.cluster.sem$cluster, c("gilman_nn_2012_cluster1.ens"="clusterI", "gilman_nn_2012_cluster2.ens"="clusterII"))

###### Mean per stage/structure ######## 
df.gilman.summary <- ddply(subset(df.gene_list, gene_list %in% gilman_lvl), c("stage", "structure_acronym"), summarise,
                                   mean = median(value, na.rm=TRUE),
                                   sd   = sd(value, na.rm=TRUE))
###### Mean per stage - FINAL ##########
df.gilman.summary.sem <- ddply(df.gilman.summary, c("stage"), summarise,
                                              mean1 = mean(mean, na.rm=TRUE),
                                              sd1   = sd(mean, na.rm=TRUE))

###################################### PROCESSING GENE lists ################################
###### Mean per stage/structure ######## 
df.summary <- ddply(df.gene_list, c("stage", "structure_acronym", "gene_list"), summarise,
                                   mean = median(value, na.rm=TRUE),
                                   sd   = sd(value, na.rm=TRUE))
## plyr magic for renaming factor level
levels(df.summary$gene_list)
df.summary$gene_list <- revalue(df.summary$gene_list, c("gene_associated.txt"="Associated Genes", "gene_nearest.txt"="Nearest Genes", "gene_prioritization.txt"="Prioritized Genes", "gene_psd_human.txt"="Post Synaptic Genes (Human)", "gene_psd_mouse.txt"="Post Synaptic Genes (Mouse)"))
levels(df.summary$gene_list)

###### Mean per stage - FINAL ##########
df.summary.sem <- ddply(df.summary, c("stage","gene_list"), summarise,
                        mean1 = mean(mean, na.rm=TRUE),
                        sd1   = sd(mean, na.rm=TRUE))


###################################### Calculating overall mean ################################
### *** Runtime ~ 10 s ***
df.all.sem <- ddply(ddply(df.expression_matrix.clean.melt, .(stage, structure_acronym), summarise, mean=median(value, na.rm=TRUE)), .(stage), summarise, mean1=mean(mean, na.rm=TRUE),  sd1=sd(mean, na.rm=TRUE))


###################################### Save DATA ################################
#save.image(file="RData_tmp/data_rnaseq_expression_processed_high_res_with_gene_lists_with_null.RData")

###################################### LOAD DATA ################################
### Whole image!
#rm(list=ls())
#load(file="RData_tmp/data_rnaseq_expression_processed_high_res_with_gene_lists_with_null.RData")

##############################################################################################################
###################################### PROCESSING GULSUNER lists - seperately ################################
###### Mean per stage/structure ######## 
str(df.gene_list)
levels(df.gene_list$gene_list)
levels(df.gene_list$structure_acronym)
#subset(df.gene_list, (gene_list=="gulsuner_S3A_damaging_cases.ens" & structure_acronym %in% c("DFC", "MFC", "OFC", "VFC")))
df.gulsuner.summary <- ddply(subset(df.gene_list, (gene_list=="gulsuner_S3A_damaging_cases.ens" & structure_acronym %in% c("DFC", "MFC", "OFC", "VFC"))), c("stage", "structure_acronym"), summarise,
                    mean = median(value, na.rm=TRUE),
                    sd   = sd(value, na.rm=TRUE))

###### Mean per stage - FINAL ##########
df.gulsuner.summary.sem <- ddply(df.gulsuner.summary, c("stage"), summarise,
                        mean1 = mean(mean, na.rm=TRUE),
                        sd1   = sd(mean, na.rm=TRUE))


###################################### PLOT GULSUNER individual genes ################################
library(hash)
file.map <- "../gene_lists/gulsuner_S3A_damaging_cases.tab"
df.map.gulsuner <- read.table(file.map, h=F)
hash.gulsuner.genes <- hash( keys=df.map.gulsuner[,1], values=df.map.gulsuner[,2] )
#df.map.gulsuner ---> col1=HGNC, col2=ENS

# Pre-natal genes
hash.gulsuner.genes[["HIF1A"]]
hash.gulsuner.genes[["CBX5"]]
# Post-natal genes
hash.gulsuner.genes[["SERPINI1"]]
# Flat/low genes
hash.gulsuner.genes[["SNX31"]]

#str(df.expression_matrix.clean.melt)
############## SINGLE PLOT - one gene
gene_hgnc <- "SNX31" # SWITCH!!
gene_ens <- as.character(hash.gulsuner.genes[[gene_hgnc]])
df.sub <- subset(df.expression_matrix.clean.melt, (ensembl_gene_id==gene_ens & structure_acronym %in% c("DFC", "MFC", "OFC", "VFC")))
p <- ggplot()
p <- p + geom_point(data=df.sub, aes(x=stage, y=value, color=structure_acronym))
p <- p + geom_smooth(data=df.sub, aes(x=stage, y=value, group=structure_acronym, color=structure_acronym), method="loess", alpha=0.1)
p + labs(title=paste(gene_hgnc,"=",gene_ens, sep=""))
#gulsuner_gene_SNX31-8x6

############## FACET wrap
#gene_hgnc <- c("HIF1A", "CBX5", "SERPINI1", "SNX31")
gene_ens <- df.map.gulsuner[,2]
df.sub <- subset(df.expression_matrix.clean.melt, (ensembl_gene_id %in% gene_ens & structure_acronym %in% c("DFC", "MFC", "OFC", "VFC")))
df.sub$hgnc <- df.map.gulsuner[match(df.sub$ensembl_gene_id, df.map.gulsuner[,2]), 1]

p <- ggplot(data=df.sub)
p <- p + geom_point(data=df.sub, aes(x=stage, y=value, color=structure_acronym))
p <- p + geom_smooth(data=df.sub, aes(x=stage, y=value, group=structure_acronym, color=structure_acronym), method="loess", alpha=0.1)
p <- p + facet_wrap(~hgnc)
p

###################################### PLOT - 1 - INCLUDING NULL ASSOCIATED ################################
########### PLOT IT! ###########

p <- ggplot()
### Adding mean Prioritized
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#d7191c', width=0.2)
p
### Adding mean ALL (df.all.sem)
p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
p
### Adding Associated Genes
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='orange', width=0.2)
#p
### Adding NULL ASSOCIATED (median)
#p <- p + geom_line(data=df.null.median.summary.sem.assoc, aes(x=stage, y=mean1, group=1, color="Associated genes (Null)"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=df.null.median.summary.sem.assoc, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='yellow', width=0.2)
#p
### Adding NULL PRIORITIZED (median)
#p <- p + geom_line(data=df.null.median.summary.sem.prio, aes(x=stage, y=mean1, group=1, color="Prioritized genes (Null)"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=df.null.median.summary.sem.prio, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='steelblue', width=0.2)
#p
### Adding GULSUNER - Specific structures "DFC", "MFC", "OFC", "VFC"
p <- p + geom_line(data=subset(df.gulsuner.summary), aes(x=stage, y=mean, group=structure_acronym, color="Gulsuner (structures)")) #linetype="Brain regions"
p
### Adding GULSUNER - mean/median over FRONTAL CORTEX ("DFC", "MFC", "OFC", "VFC")
p <- p + geom_line(data=df.gulsuner.summary.sem, aes(x=stage, y=mean1, group=1, color="Mean Gulsuner FC"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.gulsuner.summary.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
p
### Adding GULSUNER - mean/median over all structures
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "gulsuner_S3A_damaging_cases.ens"), aes(x=stage, y=mean1, group=1, color="Mean Gulsuner all structures"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "gulsuner_S3A_damaging_cases.ens"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
p



p <- p + scale_color_manual(name="Gene list", values=c("Prioritized genes (structures)"="gray", 
                                "Prioritized genes"="#d7191c",
                                "All genes"="black",
                                "Gilman et al. clusters I & II"="#2b83ba",
                                "Nearest genes"="orange",
                                "Associated genes"="orange",
                                "Post synaptic genes (mouse)"="green",
                                "Post synaptic genes (human)"="blue",
                                "Associated genes (Null)"="yellow",
                                "Prioritized genes (Null)"="steelblue",
                                guide='legend'))
p

###### Adding vertical line - prenatal vs. postnatal
p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
p

######### Adding x-tickmarks for stage
p <- do_stage_converter(p)
p

### SUPP FIG
#p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,1,0.1))))
### MAIN FIG
p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,0.1)))) #"Prioritized genes (structures)"=0.1

### VARIABLE
p <- p + labs(y="Mean brain expression")
p

