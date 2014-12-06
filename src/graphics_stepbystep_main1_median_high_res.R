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
load(file="RData/data_rnaseq_expression_processed_high_res.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt
str(df.expression_matrix.clean.melt)

######### *** NEW - 12/01/2014 ** ######## Defining new stages
# source high res stages
source("function_def_stages_high_res.R", echo=TRUE)

## sorting AGE levels
df.expression_matrix.clean.melt$age <- with(df.expression_matrix.clean.melt, factor(age, levels(age)[match(order.age, levels(age))]))
levels(df.expression_matrix.clean.melt$age)

########################################### LOAD NULL data ###################################
######### ***SKIPPED*** ########

############################# READING GENE LISTs #################################
path.datafiles <- '/Users/pascaltimshel/p_scz/brainspan/gene_lists'

###### Read into a list of files - PATTERN VERSION - read ALL .txt files in directory:
#files <- list.files(path = path.datafiles, pattern = "*.txt", full.names = TRUE) #full path
#names(files) <- list.files(path = path.datafiles, pattern = "*.txt") # filename
#cat(names(files), sep="\n")

###### Read SPECIFIC FILES:
filenames2read <- c("gene_prioritization.txt", "gene_associated.txt", "gene_nearest.txt")
filenames2read <- c(filenames2read, "gene_psd_human.txt", "gene_psd_mouse.txt")
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
#save.image(file="RData_tmp/stepbystep_median_high_res_fig1.RData")

###################################### LOAD DATA ################################
### Whole image!
#rm(list=ls())
#load(file="RData_tmp/stepbystep_median_high_res_fig1.RData")

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


###################################### STEP 0 - initialyzing ################################
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
                       "s21"="21 years",
                       "s23"="23 years",
                       "s30"="30 years",
                       "s36"="36 years",
                       "s37"="37 years",
                       "s40"="40 years")
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
  ### Prioritized (for each structure)
  p <- p + geom_line(data=subset(df.summary, gene_list == "Prioritized Genes"), aes(x=stage, y=mean, group=structure_acronym, color="Prioritized genes (structures)"), alpha=list_of_alphas[["prio_struc"]])
  ### Adding mean Prioritized
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#d7191c', width=0.2, alpha=list_of_alphas[["prio_mean"]])
  ### Adding mean ALL (df.all.sem)
  p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1, alpha=list_of_alphas[["all"]])
  p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2, alpha=list_of_alphas[["all"]])
  ### Adding Associated Genes
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='orange', width=0.2, alpha=list_of_alphas[["assoc_mean"]])
  ### Adding GILMAN to plot - MERGED cI and cII (df.gilman.summary.sem)
  p <- p + geom_line(data=df.gilman.summary.sem, aes(x=stage, y=mean1, group=1, color="Gilman et al. clusters I & II"), linetype='solid', size=1, alpha=list_of_alphas[["gilman"]]) 
  p <- p + geom_errorbar(data=df.gilman.summary.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='#2b83ba', width=0.2, alpha=list_of_alphas[["gilman"]])
  ### Adding GULSUNER - Specific structures "DFC", "MFC", "OFC", "VFC"
  p <- p + geom_line(data=subset(df.gulsuner.summary), aes(x=stage, y=mean, group=structure_acronym, color="Gulsuner et al. FC (DFC,MFC,OFC,VFC)"), alpha=list_of_alphas[["gulsuner_struc"]])
  ### Adding GULSUNER - mean/median over FRONTAL CORTEX ("DFC", "MFC", "OFC", "VFC")
  p <- p + geom_line(data=df.gulsuner.summary.sem, aes(x=stage, y=mean1, group=1, color="Gulsuner et al. FC"), linetype='solid', size=1, alpha=list_of_alphas[["gulsuner_mean_fc"]])
  p <- p + geom_errorbar(data=df.gulsuner.summary.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='cadetblue3', width=0.2, alpha=list_of_alphas[["gulsuner_mean_fc"]])
  ### Adding Nearest Genes
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, y=mean1, group=1, color="Nearest genes"), linetype='solid', size=1, alpha=list_of_alphas[["nearest"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='sky blue', width=0.2, alpha=list_of_alphas[["nearest"]])
  
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
            guide='legend')
  p <- p + scale_color_manual(name="Gene list", values=cmap)
  ##### Setting margins
  p <- p + theme(plot.margin=unit(c(5,0,0,10), "mm"))
  # ^^ margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
  #"Gulsuner et al. FC (structures)"
  #"Gulsuner et al. FC (DFC,MFC,OFC,VFC)"
    
  ### MAIN FIG
  #p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,1,0.1,1,1,0.1))))
  
  return(p)
}

#### Initialyze list of alpha
initialyze_alpha <- function() {
  names_plots <- c("prio_struc","prio_mean","all","assoc_mean","gilman","gulsuner_struc","gulsuner_mean_fc","nearest")
  n_plots <- length(names_plots)
  list_of_alphas <- as.list(rep(0,n_plots)) # vector("list", 10) ---> did not work!
  names(list_of_alphas) <- names_plots
  return(list_of_alphas)
}

#### show all plots
list_of_alphas <- lapply(initialyze_alpha(), function(x) {x <- 1})
p <- plotMe(list_of_alphas)
p

# list_of_alphas[["prio_struc"]] <- 1
# list_of_alphas[["prio_mean"]] <- 1
# list_of_alphas[["all"]] <- 1
# list_of_alphas[["assoc_mean"]] <- 1
# list_of_alphas[["gilman"]] <- 1
# list_of_alphas[["gulsuner_struc"]] <- 1
# list_of_alphas[["gulsuner_mean_fc"]] <- 1
# list_of_alphas[["nearest"]] <- 1

###################################### STEP 1 ################################
list_of_alphas <- initialyze_alpha()
### Remember to OUTCOMMENT the vertical line
p <- plotMe(list_of_alphas)
p

###################################### STEP 2 ################################
list_of_alphas <- initialyze_alpha()
### Remember to UNCOMMENT the vertical line
p <- plotMe(list_of_alphas)
p

###################################### STEP 3 - all ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,0,0,0,0,0,0), size=c(1,1,1,1,0.1,1,1,0.1))))

###################################### STEP 4 - gilman ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["gilman"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,1,0,0,0,0,0), size=c(1,1,1,1,0.1,1,1,0.1))))


###################################### STEP 5 - Gulsuner ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["gilman"]] <- 1
list_of_alphas[["gulsuner_struc"]] <- 1
list_of_alphas[["gulsuner_mean_fc"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,1,1,1,0,0,0), size=c(1,1,1,1,0.1,1,1,0.1))))


###################################### STEP 6 - associated ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["gilman"]] <- 1
list_of_alphas[["gulsuner_struc"]] <- 1
list_of_alphas[["gulsuner_mean_fc"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1

p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,1,1,1,0,0,0), size=c(1,1,1,1,0.1,1,1,0.1))))

###################################### STEP 7 - Nereast ################################

list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["gilman"]] <- 1
list_of_alphas[["gulsuner_struc"]] <- 1
list_of_alphas[["gulsuner_mean_fc"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["nearest"]] <- 1
p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,1,1,1,1,0,0), size=c(1,1,1,1,0.1,1,1,0.1))))



###################################### STEP 8 - Prioritized mean ################################

list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["gilman"]] <- 1
list_of_alphas[["gulsuner_struc"]] <- 1
list_of_alphas[["gulsuner_mean_fc"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
#list_of_alphas[["nearest"]] <- 1
list_of_alphas[["prio_mean"]] <- 1
p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,1,1,1,0,1,0), size=c(1,1,1,1,0.1,1,1,0.1))))


###################################### STEP 9 - Prioritized mean + structures ################################

list_of_alphas <- initialyze_alpha()
# list_of_alphas[["all"]] <- 1
# list_of_alphas[["gilman"]] <- 1
# list_of_alphas[["gulsuner_struc"]] <- 1
# list_of_alphas[["gulsuner_mean_fc"]] <- 1
# list_of_alphas[["assoc_mean"]] <- 1
# list_of_alphas[["nearest"]] <- 1
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_struc"]] <- 1
p <- plotMe(list_of_alphas)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(0,0,0,0,0,0,1,1), size=c(1,1,1,1,0.1,1,1,0.1))))



###################################### OPTIONAL STEPS  ################################

###################################### STEP 10 - PSD genes ################################

list_of_alphas <- initialyze_alpha()
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_struc"]] <- 1
p <- plotMe(list_of_alphas)
p

cmap <- c("Prioritized genes (structures)"="gray", 
          "Prioritized genes"="#d7191c",
          "All genes"="black",
          "Gilman et al. clusters I & II"="#2b83ba",
          "Associated genes"="darkgoldenrod1",
          "Nearest genes"="darkgoldenrod3",
          "Gulsuner et al. FC (DFC,MFC,OFC,VFC)"="cadetblue2",
          "Gulsuner et al. FC"="cadetblue3",
          "Post synaptic genes (human)"="coral1",
          "Post synaptic genes (mouse)"="coral3",
          guide='legend')

### Adding Post Synaptic Genes (Human)
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Human)"), aes(x=stage, y=mean1, group=1, color="Post synaptic genes (human)"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Human)"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='coral1', width=0.2)
### Adding Post Synaptic Genes (Mouse)
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Mouse)"), aes(x=stage, y=mean1, group=1, color="Post synaptic genes (mouse)"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Post Synaptic Genes (Mouse)"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='coral3', width=0.2)

p <- p + scale_color_manual(name="Gene list", values=cmap)
p
p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(0,0,0,0,0,0,1,1,1,1), size=c(1,1,1,1,0.1,1,1,1,1,0.1))))



