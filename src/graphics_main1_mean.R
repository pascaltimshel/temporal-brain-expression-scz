############### SYNAPSIS ###################
# Name: graphics_main1_mean.R
# This script is a newer version of "graphics_main1.R".
# The script LOADs EITHER microarray data OR PROCESSED RNAseq data.
# It was created to be able to handle multiple gene list for plotting. The script uses plyr to load in multiple gene lists.
# You may choose what files to load
# The script also process and loads GILMAN gene clusters (1, 1a, 1b and 2)
# The script was used to generate PUBLICATION READY GRAPHICS
# THIS VERSION TAKES THE MEDIAN GENE EXPRESSION VALUES.
############################################

library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)


############################# LOAD EXPRESSION DATA #################################
load(file="RData/data_marray_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
#load(file="RData/data_rnaseq_expression_processed.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt

str(df.expression_matrix.clean.melt)

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
df.gilman.summary.cluster.sem <- ddply(ddply(subset(df.gene_list, gene_list %in% gilman_lvl), .(stage, structure_acronym, gene_list), summarise, mean=mean(value, na.rm=TRUE)), .(stage, gene_list), summarise, mean1=mean(mean, na.rm=TRUE),  sd1=sd(mean, na.rm=TRUE))
## plyr magic for renaming column and factor level
df.gilman.summary.cluster.sem <- rename(df.gilman.summary.cluster.sem, c("gene_list"="cluster")) # column
df.gilman.summary.cluster.sem$cluster <- revalue(df.gilman.summary.cluster.sem$cluster, c("gilman_nn_2012_cluster1.ens"="clusterI", "gilman_nn_2012_cluster2.ens"="clusterII"))

###### Mean per stage/structure ######## 
df.gilman.summary <- ddply(subset(df.gene_list, gene_list %in% gilman_lvl), c("stage", "structure_acronym"), summarise,
                                   mean = mean(value, na.rm=TRUE),
                                   sd   = sd(value, na.rm=TRUE))
###### Mean per stage - FINAL ##########
df.gilman.summary.sem <- ddply(df.gilman.summary, c("stage"), summarise,
                                              mean1 = mean(mean, na.rm=TRUE),
                                              sd1   = sd(mean, na.rm=TRUE))

###################################### PROCESSING GENE lists ################################
###### Mean per stage/structure ######## 
df.summary <- ddply(df.gene_list, c("stage", "structure_acronym", "gene_list"), summarise,
                                   mean = mean(value, na.rm=TRUE),
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
df.all.sem <- ddply(ddply(df.expression_matrix.clean.melt, .(stage, structure_acronym), summarise, mean=mean(value, na.rm=TRUE)), .(stage), summarise, mean1=mean(mean, na.rm=TRUE),  sd1=sd(mean, na.rm=TRUE))


###################################### PLOT ################################


########### PLOT IT! ###########



p <- ggplot()
p <- p + geom_line(data=subset(df.summary, gene_list == "Prioritized Genes"), aes(x=stage, y=mean, group=structure_acronym, color="Prioritized genes (structures)")) #linetype="Brain regions"
p
### Adding mean Prioritized
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#d7191c', width=0.2)
p
### Adding mean ALL (df.all.sem)
p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
p
### Adding GILMAN to plot - MERGED cI and cII (df.gilman.summary.sem)
p <- p + geom_line(data=df.gilman.summary.sem, aes(x=stage, y=mean1, group=1, color="Gilman et al. cluster I & II"), linetype='solid', size=1) 
p <- p + geom_errorbar(data=df.gilman.summary.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='#2b83ba', width=0.2)
p
### Adding GILMAN to plot - individual cI and cII (df.gilman.summary.cluster.sem)
#p <- p + geom_line(data=df.gilman.summary.cluster.sem, aes(x=stage, y=mean1, group=cluster, linetype=cluster), color="#2b83ba", size=1) # GILMANN cluster Ia, Ib, II
#p <- p + geom_errorbar(data=df.gilman.summary.cluster.sem, aes(x=stage, group=cluster, ymax=mean1+sd1, ymin=mean1-sd1),color='#2b83ba', width=0.2)
#p
### Adding Nearest Genes
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, y=mean1, group=1, color="Nearest genes"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='sky blue', width=0.2)
#p
### Adding Associated Genes
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='orange', width=0.2)
p
### Adding Post Synaptic Genes (Human)
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Human)"), aes(x=stage, y=mean1, group=1, color="Post synaptic genes (human)"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Human)"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='black', width=0.2)
#p
### Adding Post Synaptic Genes (Mouse)
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Mouse)"), aes(x=stage, y=mean1, group=1, color="Post synaptic genes (mouse)"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Mouse)"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='black', width=0.2)
#p




p <- p + scale_color_manual(name="Gene Set", values=c("Prioritized genes (structures)"="gray", 
                                "Prioritized genes"="#d7191c",
                                "All genes"="black",
                                "Gilman et al. cluster I & II"="#2b83ba",
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
p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,1,0.1))))
### MAIN FIG
#p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,0.1)))) #"Prioritized genes (structures)"=0.1

### VARIABLE
p <- p + labs(y="Mean brain expression")
p


################### OTHER STUFF #################

### Save this for later....!
#p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(linetype=c("dashed","solid"), size=c(1.25, 0.5) )))



