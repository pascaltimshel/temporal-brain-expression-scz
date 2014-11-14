############### SYNAPSIS ###################
# This script will PLOT microarray data
# PLOT TYPES: 
# 1) Time series (including Gilman)
# 2) Correlation plot
# 3) Boxplot across structures (26)
############################################

library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

####################################### LOAD DATA ###########################################
###### LOADING EXPRESSION DATA
#source("function_read_marray.R", echo=TRUE)
# OUTPUT: df.expression_matrix.clean.melt + much more
# Loads data in "fresh" from text files
load(file="RData/data_marray_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
# gives the same is source("function_read_marray.R"), but loads data from RData

###### LOADING ADDITIONAL DATA
source("function_read_gilman.R", echo=TRUE) 
# OUTPUT: df.gilman.summary.mean



########### Setting prioritization #########
#file.gene_prioritization <- "../gene_lists/gene_prioritization.txt"
file.gene_prioritization <- "../gene_lists/gene_associated.txt"
## null - prioritized
#file.gene_prioritization <- "../gene_lists/gene_null_bmi.txt"
#file.gene_prioritization <- "../gene_lists/gene_null_height.txt"
#file.gene_prioritization <- "../gene_lists/gene_null_WHR.txt"
## null - associated
#file.gene_prioritization <- "../gene_lists/gene_null_bmi_associated.txt"
#file.gene_prioritization <- "../gene_lists/gene_null_height_associated.txt"
#file.gene_prioritization <- "../gene_lists/gene_null_WHR_associated.txt"
gene_list <- basename(file_path_sans_ext(file.gene_prioritization))

########### READ prioritization file ###########
df.gene_prioritization <- read.csv(file.gene_prioritization,h=T)
## Check mapping - compare to "df.expression_matrix.clean" NOT MOLTEN, will be too slow
sum(df.expression_matrix.clean$ensembl_gene_id %in% df.gene_prioritization[,1]) # mapped genes
sum(!df.gene_prioritization[,1] %in% df.expression_matrix.clean$ensembl_gene_id) # non-mapped genes
cat(as.character(df.gene_prioritization[!df.gene_prioritization[,1] %in% df.expression_matrix.clean$ensembl_gene_id,]), sep="\n") # print non-mapped genes

### Setting priorizied factor for MOLTEN df
df.expression_matrix.clean.melt$gene_type <- as.factor(ifelse(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.gene_prioritization[,1], "prioritized", "other"))
table(df.expression_matrix.clean.melt$gene_type)

str(df.expression_matrix.clean.melt)

################################# TIME SERIES PLOT #######################################
df.summary <- ddply(df.expression_matrix.clean.melt, c("stage", "structure_acronym", "gene_type"), summarise,
                    mean = mean(value, na.rm=TRUE),
                    sd   = sd(value, na.rm=TRUE))
str(df.summary)
df.summary.mean <- ddply(df.expression_matrix.clean.melt, c("stage", "gene_type"), summarise,
                         mean = mean(value, na.rm=TRUE),
                         sd   = sd(value, na.rm=TRUE))


p <- ggplot(subset(df.summary, gene_type == "prioritized"), aes(x=stage, y=mean, group=structure_acronym)) + geom_line(aes(colour = structure_acronym))
#p <- p + geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd), width=0.2)
p <- p + guides(colour=guide_legend(nrow = 10))
p <- p + labs(title = gene_list)
p
### Adding mean (base line)
p <- p + geom_line(data=df.summary.mean, aes(x=stage, y=mean, group = gene_type, linetype=gene_type), size=2.5)
#p <- p + guides(linetype=guide_legend(keywidth = 2, keyheight = 1))
#p <- p + geom_errorbar(data=subset(df.summary.mean, gene_type == "prioritized"), aes(x=stage, group = gene_type, linetype=gene_type, ymax = mean + sd, ymin=mean - sd), width=0.2)
p


### Adding GILMAN to plot ** GILMAN DATA MUST BE LOADED
p <- p + geom_line(data=df.gilman.summary.mean, aes(x=stage, y=mean, group = cluster, linetype=cluster), size=1.25, color="yellow")
p <- p + guides(linetype=guide_legend(keywidth = 2, keyheight = 1))
p


################################# CORRELATION PLOT #######################################

df.corr <- ddply(df.expression_matrix.clean.melt, .(stage, structure_acronym, gene_type), summarize,
                 corr=cor(gene_length, value, method="pearson", use="pairwise.complete.obs"))

df.corr.mean <- ddply(df.expression_matrix.clean.melt, .(stage, gene_type), summarize,
                      corr=cor(gene_length, value, method="pearson", use="pairwise.complete.obs"))


p <- ggplot(subset(df.corr, gene_type == "prioritized"), aes(x=stage, y=corr, group=structure_acronym)) + geom_line(aes(colour = structure_acronym))
#p <- p + geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd), width=0.2)
p <- p + guides(colour=guide_legend(nrow = 10))
p <- p + labs(title = gene_list)
p
### Adding mean (base line)
p <- p + geom_line(data=df.corr.mean, aes(x=stage, y=corr, group = gene_type, linetype=gene_type), size=2)
p <- p + guides(linetype=guide_legend(keywidth = 2, keyheight = 1))
#p <- p + geom_errorbar(data=subset(df.summary.mean, gene_type == "prioritized"), aes(x=stage, group = gene_type, linetype=gene_type, ymax = mean + sd, ymin=mean - sd), width=0.2)
p

################################# BOXPLOTs #######################################

################# TESTING using sampling
### use subset of prioritized genes
df.red <- df.expression_matrix.clean.melt[sample(1:nrow(df.expression_matrix.clean.melt), 500, replace=FALSE),]
df.red <- rbind(df.red, df.expression_matrix.clean.melt[(df.expression_matrix.clean.melt$gene_type=="prioritized"),][sample(1000,500),])

str(df.red)
p <- ggplot(df.red, aes(x=structure_acronym, y=value))
p <- p + geom_boxplot(aes(fill=gene_type)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

################# FULL
p <- ggplot(df.expression_matrix.clean.melt, aes(x=structure_acronym, y=value))
p <- p + geom_boxplot(aes(fill=gene_type)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
