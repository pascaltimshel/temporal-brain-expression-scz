############### SYNAPSIS ###################
# This script was copied from "graphics_stepbystep_main2_median_null_genes.R"
# Date of creation: 02/23/2015
# Purpose: make figures for Tune's Harvard Genetics Retreat talk
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
filenames2read <- c(filenames2read, "gene_null_bmi_associated.txt", "gene_null_bmi.txt")
filenames2read <- c(filenames2read, "gene_null_WHR_associated.txt", "gene_null_WHR.txt")

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
df.summary$gene_list <- revalue(df.summary$gene_list, c("gene_associated.txt"="SCZ Associated Genes", "gene_nearest.txt"="SCZ Nearest Genes", "gene_prioritization.txt"="SCZ Prioritized Genes"))
df.summary$gene_list <- revalue(df.summary$gene_list, c("gene_null_bmi_associated.txt"="BMI Associated Genes", "gene_null_bmi.txt"="BMI Prioritized Genes"))
df.summary$gene_list <- revalue(df.summary$gene_list, c("gene_null_WHR_associated.txt"="WHR Associated Genes", "gene_null_WHR.txt"="WHR Prioritized Genes"))
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
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "SCZ Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="SCZ Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "SCZ Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#B22600', width=0.2, alpha=list_of_alphas[["prio_mean"]])
  ### Adding Associated Genes
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "SCZ Associated Genes"), aes(x=stage, y=mean1, group=1, color="SCZ Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "SCZ Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='#FF2642', width=0.2, alpha=list_of_alphas[["assoc_mean"]])
  ### Adding mean ALL (df.all.sem)
  p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1, alpha=list_of_alphas[["all"]])
  p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2, alpha=list_of_alphas[["all"]])
  ############ NULL ############
  ### Adding NULL PRIORITIZED (median) - ribbon!
  p <- p + geom_line(data=df.null.median.summary.sem.prio, aes(x=stage, y=mean1, group=1, color="Null Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_null_line"]])
  p <- p + geom_ribbon(data=df.null.median.summary.sem.prio, aes(x=stage, group=1, ymin=mean1-sd1, ymax=mean1+sd1), alpha=list_of_alphas[["prio_null_rib"]], fill='lightskyblue3')
  p
  ### Adding NULL ASSOCIATED (median) - ribbon!
  #p <- p + geom_line(data=df.null.median.summary.sem.assoc, aes(x=stage, y=mean1, group=1, color="Null Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_null_line"]])
  #p <- p + geom_ribbon(data=df.null.median.summary.sem.assoc, aes(x=stage, group=1, ymin=mean1-sd1, ymax=mean1+sd1), alpha=list_of_alphas[["assoc_null_rib"]], fill='darkolivegreen4')
  #p

  ############ OTHER TRAITS ############
  ### BMI ###
  ### Adding BMI Prioritized
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "BMI Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="BMI Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_mean_bmi"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "BMI Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#0F2CCC', width=0.2, alpha=list_of_alphas[["prio_mean_bmi"]])
  ### Adding BMI Associated
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "BMI Associated Genes"), aes(x=stage, y=mean1, group=1, color="BMI Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_mean_bmi"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "BMI Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='#66C1FF', width=0.2, alpha=list_of_alphas[["assoc_mean_bmi"]])
  
  ### WHR ###
  ### Adding WHR Prioritized
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "WHR Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="WHR Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_mean_whr"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "WHR Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#5FB21F', width=0.2, alpha=list_of_alphas[["prio_mean_whr"]])
  ### Adding WHR Associated
  #p <- p + geom_line(data=subset(df.summary.sem, gene_list == "WHR Associated Genes"), aes(x=stage, y=mean1, group=1, color="WHR Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_mean_whr"]])
  #p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "WHR Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='#96FF46', width=0.2, alpha=list_of_alphas[["assoc_mean_whr"]])
  
  
  
  ###### Adding vertical line - prenatal vs. postnatal
  #p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
  
  ######### Adding title + x-tickmarks
  p <- do_stage_converter(p)
  p
  #### SETTING LEGEND
  cmap <- c("SCZ Prioritized genes"="#B22600", # #d7191c
            "SCZ Associated genes"="#FF2642", #"#FF622C",
            "All genes"="black",
            "Null Associated genes"="darkolivegreen4",
            "Null Prioritized genes"="lightskyblue3",
            "BMI Prioritized genes"="#0F2CCC",
            "BMI Associated genes"= "#66C1FF", #"#4B5BB2",
            "WHR Prioritized genes"="#5FB21F",
            "WHR Associated genes"="#96FF46",
            guide='legend')
  
#   cmap <- c("SCZ Prioritized genes"="#B22600", # #d7191c
#             "SCZ Associated genes"="darkgoldenrod1",
#             "All genes"="black",
#             "Null Associated genes"="darkolivegreen4",
#             "Null Prioritized genes"="lightskyblue3",
#             "BMI Prioritized genes"="cadetblue2",
#             "BMI Associated genes"="cadetblue3",
#             "WHR Prioritized genes"="black",
#             "WHR Associated genes"="black",
#             guide='legend')
  
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
  names_plots <- c(names_plots, "prio_mean_bmi", "assoc_mean_bmi") # NEW 2015: BMI
  names_plots <- c(names_plots, "prio_mean_whr", "assoc_mean_whr") # NEW 2015: WHR
  n_plots <- length(names_plots)
  list_of_alphas <- as.list(rep(0,n_plots)) # vector("list", 10) ---> did not work!
  names(list_of_alphas) <- names_plots
  return(list_of_alphas)
}


#### show all plots
list_of_alphas <- lapply(initialyze_alpha(), function(x) {x <- 1})
list_of_alphas
p <- plotMe(list_of_alphas)
p

# list_of_alphas[["prio_mean"]] <- 1
# list_of_alphas[["all"]] <- 1
# list_of_alphas[["assoc_mean"]] <- 1
# list_of_alphas[["prio_null_line"]] <- 1
# list_of_alphas[["prio_null_rib"]] <- 1
# list_of_alphas[["assoc_null_line"]] <- 1
# list_of_alphas[["assoc_null_rib"]] <- 1

###################################### PLOT PARAMS ################################
plot.param.width = 10
plot.param.height = 6

###################################### STEP 0a - blank ################################
list_of_alphas <- initialyze_alpha()
p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=0)))
ggsave("step-0a_blank.pdf", w=plot.param.width, h=plot.param.height)

###################################### STEP 0b - baseline (all) ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,0,0,0,0,0))))
ggsave("step-0b_baseline.pdf", w=plot.param.width, h=plot.param.height)

###################################### STEP 1a - assoc_mean [SCZ] ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,0,0,1,0,0))))
ggsave("step-1a_assoc.pdf", w=plot.param.width, h=plot.param.height)

###################################### STEP 1a - assoc_mean [SCZ + BMI] ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_mean_bmi"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0,0,1,0,0))))
ggsave("step-1b_assoc.pdf", w=plot.param.width, h=plot.param.height)

###################################### STEP 2a - assoc_mean [SCZ + BMI] + prio_mean [SCZ] ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_mean_bmi"]] <- 1
list_of_alphas[["prio_mean"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0,0,1,1,0))))
ggsave("step-2a_assoc_prio.pdf", w=plot.param.width, h=plot.param.height)

###################################### STEP 2b - assoc_mean [SCZ + BMI] + prio_mean [SCZ + BMI] ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_mean_bmi"]] <- 1
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_mean_bmi"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,1,0,1,1,0))))
ggsave("step-2b_assoc_prio.pdf", w=plot.param.width, h=plot.param.height)

###################################### STEP 4 - prio_mean [SCZ + BMI + WHR] + NULL ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
#list_of_alphas[["assoc_mean"]] <- 1
#list_of_alphas[["assoc_mean_bmi"]] <- 1
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_mean_bmi"]] <- 1
list_of_alphas[["prio_mean_whr"]] <- 1

list_of_alphas[["prio_null_line"]] <- 1
list_of_alphas[["prio_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,1,0.1,0,1,1), size=c(1,1,1,4,1,1,1))))
ggsave("step-3a_prio_null.pdf", w=plot.param.width, h=plot.param.height)












