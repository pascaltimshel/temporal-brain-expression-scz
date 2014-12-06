library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)



########################################### LOAD data ###################################

############################# EXTRACT BROAD DATA #################################
source("function_GWAScatalog_get_summary.R")
################################ ** PRIORITIZED genes ** #########################
load("RData/GWAScatalog_prioritized_RData_broad_rnaseq_priority.RData")
names(list.par_analysis)
############### Extracting from list
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal.paired <- llply(list.null.fits, "[[", "fit.natal.paired") # *NEW PAIRED* 


############################# SUBSETTING on criteria #################################
####### csv file to filter on the number of genes - ONLY THE PRIORITIZED GENES SHOULD BE USED FOR SUBSETTING.
file.gwas_catalog_stat <- "../data/GWAS_catalog_prioritized_genes.csv"
#file.gwas_catalog_stat <- "../data/GWAS_catalog_associated_genes.csv"

####### read csv file
df.gwas_catalog_stat <- read.csv(file.gwas_catalog_stat)
str(df.gwas_catalog_stat)
df.gwas_catalog_stat.summary <- ddply(df.gwas_catalog_stat, .(permutation), summarise,
                                      n_genes = length(n_genes))
criteria <- 10
df.gwas_catalog_stat.summary.crit <- subset(df.gwas_catalog_stat.summary, n_genes>=criteria)
phenotypes_passing_criteria <- as.character(df.gwas_catalog_stat.summary.crit$permutation)

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

############################# FACET PLOT: plot #################################
### Load image from "stepbystep_median_null.R" in "src_presentation/"
#rm(list=ls())
#load(file="RData/stepbystep_median_null_RNASEQ.RData") # --> image file
#save(df.summary, df.summary.sem, df.all.sem, file="RData/graphics_GWAS-catalog_facet_plot.RData")
load(file="RData/graphics_GWAS-catalog_facet_plot.RData")

########### PREPARING SUBSET ##########
### subsetting df.null.median.*
# prio
df.null.median.summary.prio.sub <- subset(df.null.median.summary.prio, permutation %in% phenotypes_passing_criteria)
df.null.median.summary.sem.prio.sub <- subset(df.null.median.summary.sem.prio, permutation %in% phenotypes_passing_criteria)
# assoc
df.null.median.summary.assoc.sub <- subset(df.null.median.summary.assoc, permutation %in% phenotypes_passing_criteria)
df.null.median.summary.sem.assoc.sub <- subset(df.null.median.summary.sem.assoc, permutation %in% phenotypes_passing_criteria)

####### MODIFYING phenotype names #######
## IMPORTANT: *substitute underscores ("_") in phenotype name*##
df.null.median.summary.prio.sub$permutation <- gsub("_", " ", df.null.median.summary.prio.sub$permutation)
df.null.median.summary.sem.prio.sub$permutation <- gsub("_", " ", df.null.median.summary.sem.prio.sub$permutation)
df.null.median.summary.assoc.sub$permutation <- gsub("_", " ", df.null.median.summary.assoc.sub$permutation)
df.null.median.summary.sem.assoc.sub$permutation <- gsub("_", " ", df.null.median.summary.sem.assoc.sub$permutation)

## adding "wrapped" text column ##
#Metabolic_traits --> n_characters=16
wrapit <- function(x) {paste(strwrap(x,width=26),collapse="\n")}
# ^^ width=26 works for the following parameters: strip.text.x=element_text(size=16), save as pdf 20x20
df.null.median.summary.prio.sub$wrap_permutation <- sapply(df.null.median.summary.prio.sub$permutation, wrapit)
df.null.median.summary.sem.prio.sub$wrap_permutation <- sapply(df.null.median.summary.sem.prio.sub$permutation, wrapit)
df.null.median.summary.assoc.sub$wrap_permutation <- sapply(df.null.median.summary.assoc.sub$permutation, wrapit)
df.null.median.summary.sem.assoc.sub$wrap_permutation <- sapply(df.null.median.summary.sem.assoc.sub$permutation, wrapit)



str(df.summary.sem)
str(df.null.median.summary.prio.sub)
str(df.null.median.summary.sem.prio.sub)

########### PLOT IT ##########
p <- ggplot(data=df.null.median.summary.sem.prio.sub)
### SCZ - priotitized genes
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="SCZ prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='gray', width=0.2)
### GWAScatalog - prioritized (for each structure)
p <- p + geom_line(data=df.null.median.summary.prio.sub, aes(x=stage, y=mean, group=structure_acronym, color="GWAS catalog prioritized genes (structures)"))
### GWAScatalog - mean prioritized
p <- p + geom_line(data=df.null.median.summary.sem.prio.sub, aes(x=stage, y=mean1, group=1, color="GWAS catalog prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.null.median.summary.sem.prio.sub, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
### GWAScatalog - mean associated
p <- p + geom_line(data=df.null.median.summary.sem.assoc.sub, aes(x=stage, y=mean1, group=1, color="GWAS catalog associated genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.null.median.summary.sem.assoc.sub, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
#p <- p + facet_wrap(~permutation)
p <- p + facet_wrap(~wrap_permutation)
#p

cmap <- c("SCZ prioritized genes"="gray", 
          "GWAS catalog prioritized genes (structures)"="chartreuse3",
          "GWAS catalog prioritized genes"="brown1",
          "GWAS catalog associated genes"="cyan3",
          guide='legend')
p <- p + scale_color_manual(name="Gene list", values=cmap)
p <- do_stage_converter(p)
## Setting x-axis text size (DIFFICULT!)
p <- p + theme(axis.text.x=element_text(size=8))
p <- p + labs(y="Mean brain expression")
### Setting text size for label/strip.text in facets
p <- p + theme(strip.text.x=element_text(size=16, colour="black"))
### Setting axis.title size for y-axis
p <- p + theme(axis.title.y=element_text(size=16))
p

#GWAS-catalog_facet_min10_prioritized_associated_incl_SCZ-20x20.pdf
#facet_min10_prioritized_associated_incl_SCZ-25x25.pdf

