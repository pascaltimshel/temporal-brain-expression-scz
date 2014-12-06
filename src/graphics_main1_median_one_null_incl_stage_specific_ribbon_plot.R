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


############################# LOAD EXPRESSION DATA #################################
#load(file="RData/data_marray_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
load(file="RData/data_rnaseq_expression_processed.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt

str(df.expression_matrix.clean.melt)

########################################### LOAD NULL data ###################################
############ ** ASSOCIATED genes ** ##########
### marray
#load("RData/null_RData_broad_marray_associated_priority.RData") #time_elapsed, list.par_analysis
### rnaseq
#load("RData/null_RData_broad_rnaseq_associated_priority.RData") #time_elapsed, list.par_analysis

############ ** PRIORITIZED genes ** ##########
### marray
#load("RData/null_RData_broad_marray_prioritized_priority.RData") #time_elapsed, list.par_analysis
### rnaseq
load("RData/null_RData_broad_rnaseq_prioritized_priority.RData") #time_elapsed, list.par_analysis


############################# EXTRACT BROAD DATA #################################
############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal <- lapply(list.null.fits, "[[", "fit.natal")
list.null.fit.prioritized.higher <- lapply(list.null.fits, "[[", "fit.prioritized.higher")
list.null.list.fit.stage <- lapply(list.null.fits, "[[", "list.fit.stage")
############## Extracting from list.null.list.fit.stage
df.fit.stage<-ldply(seq_along(list.null.list.fit.stage), function(i) {ldply(list.null.list.fit.stage[[i]], function(fit.stage) {data.frame(t_statistic=fit.stage$statistic, perm=i)})})

############## Combining summary data frames
df.null.mean.summary <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary <- ldply(list.null.median.summary) # COMBINING list of data frames
df.null.median.summary.per_perm <- ddply(df.null.median.summary, .(permutation, stage), summarise,
                                              mean_stage_across_structure = mean(mean, na.rm=TRUE),
                                              sd_stage_across_structure   = sd(mean, na.rm=TRUE))
df.null.median.summary.sem <- ddply(df.null.median.summary.per_perm, .(stage), summarise,
                                         mean1 = mean(mean_stage_across_structure, na.rm=TRUE),
                                         sd1   = sd(mean_stage_across_structure, na.rm=TRUE))

############################# READING GENE LISTs #################################
path.datafiles <- '/Users/pascaltimshel/p_scz/brainspan/gene_lists'

###### Read into a list of files - PATTERN VERSION - read ALL .txt files in directory:
#files <- list.files(path = path.datafiles, pattern = "*.txt", full.names = TRUE) #full path
#names(files) <- list.files(path = path.datafiles, pattern = "*.txt") # filename
#cat(names(files), sep="\n")

###### Read SPECIFIC FILES:
filenames2read <- c("gene_associated.txt", "gene_nearest.txt", "gene_prioritization.txt", "gene_psd_human.txt", "gene_psd_mouse.txt")
filenames2read <- c(filenames2read, "gilman_nn_2012_cluster1.ens", "gilman_nn_2012_cluster1a.ens", "gilman_nn_2012_cluster1b.ens", "gilman_nn_2012_cluster2.ens")
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
#save(file="RData/data_rnaseq_expression_processed_with_gene_lists.RData")





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
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='orange', width=0.2)
p
### Adding NULL ASSOCIATED (median)
#qplot(stage, mean1, group=permutation, data=df.null.median.summary.sem, geom="line")
p <- p + geom_line(data=df.null.median.summary.sem, aes(x=stage, y=mean1, group=1, color="Associated genes (Null)"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.null.median.summary.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='orange', width=0.2)
p + labs(title="???") #labs(title="RNAseq")
#rnaseq_incl_associated_null-9x6

p <- p + scale_color_manual(name="Gene list", values=c("Prioritized genes (structures)"="gray", 
                                "Prioritized genes"="#d7191c",
                                "All genes"="black",
                                "Gilman et al. clusters I & II"="#2b83ba",
                                "Nearest genes"="orange",
                                "Associated genes"="orange",
                                "Post synaptic genes (mouse)"="green",
                                "Post synaptic genes (human)"="blue",
                                guide='legend'))
p

###### Adding vertical line - prenatal vs. postnatal
p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
p


######### Adding x-tickmarks for stage
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
p


### SUPP FIG
#p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,1,0.1))))
### MAIN FIG
p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,0.1)))) #"Prioritized genes (structures)"=0.1

### VARIABLE
p <- p + labs(y="Mean brain expression")
p



###################################### PLOT - 2 ################################
# *** MUST RUN/LOAD "T-test: test for each specific time point" for the CORRECT DATA SET (RNAseq/Microarray) before performning the following
#### LOADING: df.fit.stage.obs
### RNAseq
load(file="RData/fit_stage_obs_rnaseq_associated.RData")
df.fit.stage.obs.rnaseq.associated <- df.fit.stage.obs
load(file="RData/fit_stage_obs_rnaseq_prioritized.RData")
df.fit.stage.obs.rnaseq.prioritized <- df.fit.stage.obs
### Microarray
load(file="RData/fit_stage_obs_marray_prioritized.RData")
df.fit.stage.obs.marray.prioritized <- df.fit.stage.obs

# df.fit.stage.obs
p <- ggplot()
### Adding observed t-statistic df.fit.stage.obs - ASSOCIATED
p <- p + geom_line(data=df.fit.stage.obs.rnaseq.associated, aes(x=stage, y=t_statistic, group=1, color="Associated"), size=1)
p <- p + geom_line(data=df.fit.stage.obs.rnaseq.prioritized, aes(x=stage, y=t_statistic, group=1, color="Prioritized"), size=1)
p <- p + geom_line(data=df.fit.stage.obs.marray.prioritized, aes(x=stage, y=t_statistic, group=1, color="Microarray - Prioritized"), size=1)
p
### plot loess smooting
#p + stat_smooth(data=df.fit.stage, aes(x=stage, y=t_statistic, group=1), method="loess",se=T, level=0.95)
### plot sd as ribbon
df.fit.stage.summary <- ddply(df.fit.stage, .(stage), summarise, mean=mean(t_statistic), sd=sd(t_statistic))
p <- p + geom_line(data=df.fit.stage.summary, aes(x=stage, y=mean, group=1, color="Null genes (associated, RNAseq)"), size=1)
p <- p + geom_ribbon(data=df.fit.stage.summary, aes(x=stage, group = 1, ymin=mean-sd, ymax=mean+sd), alpha=0.5, fill="black", color="gray")
p + labs(title="RNAseq - ASSOCIATED GENES - observed and null t-statistic")



