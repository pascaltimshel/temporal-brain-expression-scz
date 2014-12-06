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


############################# LOAD EXPRESSION DATA #################################
#load(file="RData/data_marray_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
load(file="RData/data_rnaseq_expression_processed.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt
str(df.expression_matrix.clean.melt)

########################################### LOAD NULL data ###################################
################################## ** ASSOCIATED genes ** ##########################
### marray
#load("RData/null_RData_broad_marray_associated_priority.RData") #time_elapsed, list.par_analysis
### rnaseq
load("RData/null_RData_broad_rnaseq_associated_priority.RData") #time_elapsed, list.par_analysis

############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
############## Combining summary data frames
df.null.mean.summary.assoc <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary.assoc <- ldply(list.null.median.summary) # COMBINING list of data frames
# THIS BELOW APROACH WAS USED TO GENERATE PLOT FOR *HIRCHHORN LAB PRESENTATION*
# df.null.median.summary.sem.assoc <- ddply(df.null.median.summary.assoc, c("stage"), summarise,
#                                          mean1 = mean(mean, na.rm=TRUE),
#                                          sd1   = sd(mean, na.rm=TRUE))
####
# First we find the stage average across structures FOR EACH PERMUTATION. 
# In this way we process each permutation as we did for SCZ prioritized genes
df.null.median.summary.per_perm.assoc <- ddply(df.null.median.summary.assoc, .(permutation, stage), summarise,
                                          mean_stage_across_structure = mean(mean, na.rm=TRUE),
                                          sd_stage_across_structure   = sd(mean, na.rm=TRUE))

df.null.median.summary.sem.assoc <- ddply(df.null.median.summary.per_perm.assoc, .(stage), summarise,
                                          mean1 = mean(mean_stage_across_structure, na.rm=TRUE),
                                          sd1   = sd(mean_stage_across_structure, na.rm=TRUE))



################################ ** PRIORITIZED genes ** #########################
### marray
#load("RData/null_RData_broad_marray_prioritized_priority.RData") #time_elapsed, list.par_analysis
### rnaseq
load("RData/null_RData_broad_rnaseq_prioritized_priority.RData") #time_elapsed, list.par_analysis

############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
############## Combining summary data frames
df.null.mean.summary.prio <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary.prio <- ldply(list.null.median.summary) # COMBINING list of data frames
df.null.median.summary.per_perm.prio <- ddply(df.null.median.summary.prio, .(permutation, stage), summarise,
                                               mean_stage_across_structure = mean(mean, na.rm=TRUE),
                                               sd_stage_across_structure   = sd(mean, na.rm=TRUE))

df.null.median.summary.sem.prio <- ddply(df.null.median.summary.per_perm.prio, .(stage), summarise,
                                          mean1 = mean(mean_stage_across_structure, na.rm=TRUE),
                                          sd1   = sd(mean_stage_across_structure, na.rm=TRUE))


############################# READING GENE LISTs #################################
path.datafiles <- '/Users/pascaltimshel/p_scz/brainspan/gene_lists'

###### Read into a list of files - PATTERN VERSION - read ALL .txt files in directory:
#files <- list.files(path = path.datafiles, pattern = "*.txt", full.names = TRUE) #full path
#names(files) <- list.files(path = path.datafiles, pattern = "*.txt") # filename
#cat(names(files), sep="\n")

###### Read SPECIFIC FILES:
filenames2read <- c("gene_prioritization.txt", "gene_associated.txt", "gene_nearest.txt")
#filenames2read <- c(filenames2read, "gene_psd_human.txt", "gene_psd_mouse.txt")
#filenames2read <- c(filenames2read, "gilman_nn_2012_cluster1.ens", "gilman_nn_2012_cluster1a.ens", "gilman_nn_2012_cluster1b.ens", "gilman_nn_2012_cluster2.ens")
#filenames2read <- c(filenames2read, "gulsuner_S3A_damaging_cases.ens")
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
# *** REMOVED ** #

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
#save.image(file="RData/stepbystep_median_null_RNASEQ.RData")
#save.image(file="RData/stepbystep_median_null_MARRAY.RData")

###################################### LOAD DATA ################################
#rm(list=ls())
#load(file="RData/stepbystep_median_null_RNASEQ.RData") # --> image file


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
  p <- p + scale_x_discrete(name="", labels = stage_converter) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  p <- p + theme(axis.title.y = element_text(size=15)) # axis title (e.g. "Mean brain expression")
  p <- p + theme(axis.text.y = element_text(size=15))
  p <- p + theme(axis.text.x = element_text(size=12))
  p <- p + labs(y="Mean brain expression")
  return(p)
}


###################################### STEP 0 - final plot ################################
plotMe <- function(list_of_alphas) {
  library("grid") # needed for unit() function
  
  p <- ggplot()
  ### Adding mean Prioritized
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#d7191c', width=0.2, alpha=list_of_alphas[["prio_mean"]])
  ### Adding mean ALL (df.all.sem)
  p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1, alpha=list_of_alphas[["all"]])
  p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2, alpha=list_of_alphas[["all"]])
  ### Adding Associated Genes
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='darkgoldenrod1', width=0.2, alpha=list_of_alphas[["assoc_mean"]])
  ### Adding NULL ASSOCIATED (median) - ribbon!
  p <- p + geom_line(data=df.null.median.summary.sem.assoc, aes(x=stage, y=mean1, group=1, color="Associated genes (Null)"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_null_line"]])
  p <- p + geom_ribbon(data=df.null.median.summary.sem.assoc, aes(x=stage, group=1, ymin=mean1-sd1, ymax=mean1+sd1), alpha=list_of_alphas[["assoc_null_rib"]], fill='darkolivegreen4')
  p
  ### Adding NULL PRIORITIZED (median) - ribbon!
  p <- p + geom_line(data=df.null.median.summary.sem.prio, aes(x=stage, y=mean1, group=1, color="Prioritized genes (Null)"), linetype='solid', size=1, alpha=list_of_alphas[["prio_null_line"]])
  p <- p + geom_ribbon(data=df.null.median.summary.sem.prio, aes(x=stage, group=1, ymin=mean1-sd1, ymax=mean1+sd1), alpha=list_of_alphas[["prio_null_rib"]], fill='lightskyblue3')
  p

  ###### Adding vertical line - prenatal vs. postnatal
  p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
  
  ######### Adding title + x-tickmarks
  p <- do_stage_converter(p)
  p
  #### SETTING LEGEND
  cmap <- c("Prioritized genes (structures)"="gray", 
            "Prioritized genes"="#d7191c",
            "All genes"="black",
            "Gilman et al. clusters I & II"="#2b83ba",
            "Associated genes"="darkgoldenrod1",
            "Nearest genes"="darkgoldenrod3",
            "Gulsuner et al. FC (DFC,MFC,OFC,VFC)"="cadetblue2",
            "Gulsuner et al. FC"="cadetblue3",
            "Associated genes (Null)"="darkolivegreen4",
            "Prioritized genes (Null)"="lightskyblue3",
            guide='legend')
  
  p <- p + scale_color_manual(name="Gene list", values=cmap)
  ##### Setting margins
  p <- p + theme(plot.margin=unit(c(5,0,0,10), "mm"))
  # ^^ margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
  

  

  ### MAIN FIG
  #p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,1,0.1,1,1,0.1))))
  
  return(p)
}


#### Initialyze list of alpha
initialyze_alpha <- function() {
  names_plots <- c("prio_mean","all","assoc_mean")
  names_plots <- c(names_plots, "prio_null_line","prio_null_rib","assoc_null_line","assoc_null_rib")
  n_plots <- length(names_plots)
  list_of_alphas <- as.list(rep(0,n_plots)) # vector("list", 10) ---> did not work!
  names(list_of_alphas) <- names_plots
  return(list_of_alphas)
}


#### show all plots
list_of_alphas <- lapply(initialyze_alpha(), function(x) {x <- 1})
p <- plotMe(list_of_alphas)
p

# list_of_alphas[["prio_mean"]] <- 1
# list_of_alphas[["all"]] <- 1
# list_of_alphas[["assoc_mean"]] <- 1
# list_of_alphas[["prio_null_line"]] <- 1
# list_of_alphas[["prio_null_rib"]] <- 1
# list_of_alphas[["assoc_null_line"]] <- 1
# list_of_alphas[["assoc_null_rib"]] <- 1


###################################### STEP 1 - assoc_mean ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,0,0,0), size=c(1,1,3,1,3))))
#p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0.1,1,0.1), size=c(1,1,3,1,3))))

###################################### STEP 2 - assoc_mean + null ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_null_line"]] <- 1
list_of_alphas[["assoc_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0.1,0,0), size=c(1,1,3,1,3))))

###################################### STEP 3 - prio_mean + null ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_null_line"]] <- 1
list_of_alphas[["prio_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,0,1,0.1), size=c(1,1,3,1,3))))

###################################### STEP 4 - all ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_null_line"]] <- 1
list_of_alphas[["assoc_null_rib"]] <- 0.3
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_null_line"]] <- 1
list_of_alphas[["prio_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0.1,1,0.1), size=c(1,1,3,1,3))))


###################################### OPTIONAL ################################
###################################### STEP 5 - MICROARRAY ################################
#rm(list=ls())
#load(file="RData/stepbystep_median_null_MARRAY.RData") # --> image file
################ LOAD DATA ##########

list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_null_line"]] <- 1
list_of_alphas[["assoc_null_rib"]] <- 0.3
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_null_line"]] <- 1
list_of_alphas[["prio_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p <- p + labs(y="Mean brain expression", title="Microarray")
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0.1,1,0.1), size=c(1,1,3,1,3))))
### VARIABLE

















