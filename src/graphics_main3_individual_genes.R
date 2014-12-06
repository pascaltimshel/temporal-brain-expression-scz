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
library(grid) # for unit() function


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)


############################# LOAD EXPRESSION DATA - HIGH RES #################################
#load(file="RData/data_marray_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
load(file="RData/data_rnaseq_expression_processed_full_res.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt
str(df.expression_matrix.clean.melt)

######## Defining new stages #########
source("function_def_stages.R", echo=TRUE)
source("function_load_gene_list_data_individual_gene_plots.R", echo=TRUE)
  # --> do_stage_converter
  # --> df.summary
  # --> df.summary.sem
  # --> df.all.sem

###################################### Save DATA ################################
#save.image(file="RData_tmp/XXXXXXXXXXX.RData")

########## Read ENSEMBL-2-HGNC mapping file for prioritized genes ##############
file.map.prioritized <- '/Users/pascaltimshel/p_scz/brainspan/gene_lists/gene_prioritization.csv'
df.map.prioritized <- read.csv(file.map.prioritized, h=F)

###################################### PLOT individual genes ################################

############## SINGLE PLOT - one gene
### Subsetting genes
#gene_ens <- as.character(list_of_data[["gene_prioritization.txt"]][1,]) # OLD METHOD: ONE GENE
gene_ens <- as.character(df.map.prioritized[1,1]) # ONE GENE
df.sub <- subset(df.expression_matrix.clean.melt, (ensembl_gene_id==gene_ens & structure_acronym %in% c("DFC", "MFC", "OFC", "VFC")))
### constructing label: mapping ensembl-->hgnc
row_idx <- match(df.sub$ensembl_gene_id, df.map.prioritized[,1]) # No need for as.character(): match() --> factors are converted to character vectors,
hgnc_symbol <- as.character(df.map.prioritized[row_idx,2])
replacement_string <- paste(df.sub$ensembl_gene_id, "=", hgnc_symbol, sep="")
df.sub$facet_label <- replacement_string
### Making plot
p <- ggplot()
p <- p + geom_point(data=df.sub, aes(x=stage, y=value, color=structure_acronym))
p <- p + geom_smooth(data=df.sub, aes(x=stage, y=value, group=structure_acronym, color=structure_acronym), method="loess", alpha=0.1)
p + labs(title=unique(df.sub$facet_label))


############## FACET wrap
### Subsetting genes
gene_ens <- as.character(df.map.prioritized[,1]) # ALL GENES!
df.sub <- subset(df.expression_matrix.clean.melt, (ensembl_gene_id %in% gene_ens & structure_acronym %in% c("DFC", "MFC", "OFC", "VFC")))
### constructing label: mapping ensembl-->hgnc
row_idx <- match(df.sub$ensembl_gene_id, df.map.prioritized[,1]) # No need for as.character(): match() --> factors are converted to character vectors,
hgnc_symbol <- as.character(df.map.prioritized[row_idx,2])
# replacing unmapped genes (hgnc=="") with their ENSEMBL names:
df.sub$facet_label <- ifelse(hgnc_symbol!="", hgnc_symbol, as.character(df.sub$ensembl_gene_id))
### Making plot
p <- ggplot(data=df.sub)
p <- p + geom_point(data=df.sub, aes(x=stage, y=value, color=structure_acronym))
p <- p + geom_smooth(data=df.sub, aes(x=stage, y=value, group=structure_acronym, color=structure_acronym), method="loess", alpha=0.1)
p <- p + facet_wrap(~ facet_label)
## add mean of prioritized
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Mean all structures"), linetype='solid', size=1)
## add mean of all genes
#p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1), linetype='solid', color="Black", size=1)
#### SETTING LEGEND
cmap <- c("Mean all structures"="gray",
          "DFC"="#F8766D",
          "MFC"="#7CAE00",
          "OFC"="#00BFC4",
          "VFC"="#C77CFF",
          guide='legend')
p <- p + scale_color_manual(name="Brain structure", values=cmap)
p <- do_stage_converter(p)
## Setting x-axis text size (DIFFICULT!)
p <- p + theme(axis.text.x=element_text(size=8))
p <- p + labs(y="Mean brain expression")
### Setting text size for label/strip.text in facets
p <- p + theme(strip.text.x=element_text(size=16, colour="black", face="italic"))
### Setting axis.title size for y-axis
p <- p + theme(axis.title.y=element_text(size=20))
### Setting longer distance between x-axis tick marks
#p <- p + scale_x_discrete(expand=c(0.2,0)) # DOES NOT WORK. multiplicative and additive expansion constants
#p <- p + theme(axis.ticks.length=unit(10, 'mm')) #DOES NOT WORK.
#p + scale_x_discrete(breaks=c(1,2,3), labels=c("A1","A2","A3")) #DOES NOT WORK.
p
#prioritized_genes_structure=FC_with_mean-20x20.pdf



################################### HELPER FUNCTIONS ###############################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

cols = gg_color_hue(4)
cols
scales:::show_col(cols)


########## PREVIOUS COLORS ##########
# 4-class RdBu
# http://colorbrewer2.org/?type=diverging&scheme=RdBu&n=4

# cmap <- c("Mean all structures"="gray",
#           "DFC"="#d7191c",
#           "MFC"="#fdae61",
#           "OFC"="#abd9e9",
#           "VFC"="#2c7bb6",
#           guide='legend')
