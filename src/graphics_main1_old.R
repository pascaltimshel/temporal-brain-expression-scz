############### SYNAPSIS ###################
# This script is the "old" main graphic generator.
# The script loads either Microarray (marray) or RNAseq data
# IMPORTANTLY: this script uses the old "setting prioritization" scheme
# i.e it cannot handle multiple gene list (e.g. associated and prioritized genes) at the same same.
# IMPORTANTLY: this script contains "leftovers" of ggplot LEGEND experimentation

# DEPENDENCIES: this script is dependent on loading the Gilman RData file (df.gilman.melt)
############################################

library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)


############################# LOAD EXPRESSION DATA #################################
load("data_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
load("data_Gilman.RData") # df.gilman.melt | levels(df.gilman.melt$cluster) --> "cluster1a" "cluster1b" "cluster2"


############ RNAseq
#load(file="data_rnaseq_expression.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt
#### USE THIS - after processing:
#load("data_rnaseq_Gilman.RData") # RNAseq AFTER PROCESSING | df.gilman.melt | levels(df.gilman.melt$cluster) --> "cluster1a" "cluster1b" "cluster2"
#load(file="data_rnaseq_expression_processed.RData") # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt


########### Setting prioritization #########
file.gene_prioritization <- "/Users/pascaltimshel/p_scz/brainspan/data/141031/microarray/gene_prioritization.txt"
#file.gene_prioritization <- "/Users/pascaltimshel/p_scz/brainspan/data/141031/microarray/gene_associated.txt"
#file.gene_prioritization <- "/Users/pascaltimshel/p_scz/brainspan/data/141031/microarray/gene_nearest.txt"
#file.gene_prioritization <- "/Users/pascaltimshel/p_scz/brainspan/data/141031/microarray/gene_PSD.txt"


gene_list <- basename(file_path_sans_ext(file.gene_prioritization))

########### READ prioritization file ###########
df.gene_prioritization <- read.csv(file.gene_prioritization,h=T)

### Setting priorizied factor for MOLTEN df
df.expression_matrix.clean.melt$gene_type <- as.factor(ifelse(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.gene_prioritization[,1], "prioritized", "other"))
table(df.expression_matrix.clean.melt$gene_type)
str(df.expression_matrix.clean.melt)

############ TIME SERIES PLOT ##############
df.summary <- ddply(df.expression_matrix.clean.melt, c("stage", "structure_acronym", "gene_type"), summarise,
                    mean = mean(value, na.rm=TRUE),
                    sd   = sd(value, na.rm=TRUE))
str(df.summary)
### STANDARD ERROR OF THE MEAN - *OBS name change! mean1, sd1
df.summary.mean <- ddply(df.summary, c("stage", "gene_type"), summarise,
                         mean1 = mean(mean, na.rm=TRUE),
                         sd1   = sd(mean, na.rm=TRUE)) # this is S.E.M
str(df.summary.mean)

########## GILMAN data manipulation
str(df.gilman.melt)
levels(df.gilman.melt$cluster)
df.gilman.summary <- ddply(df.gilman.melt, c("stage", "structure_acronym"), summarise,
                                mean = mean(value, na.rm=TRUE),
                                sd   = sd(value, na.rm=TRUE))

df.gilman.summary.cluster <- ddply(df.gilman.melt, c("stage", "cluster"), summarise,
                           mean = mean(value, na.rm=TRUE),
                           sd   = sd(value, na.rm=TRUE))

df.gilman.summary.mean.clustermerged <- ddply(df.gilman.summary, c("stage"), summarise,
                                mean1 = mean(mean, na.rm=TRUE),
                                sd1   = sd(mean, na.rm=TRUE))

########### PLOT IT! ###########

p <- ggplot() 
p <- p + geom_line(data=subset(df.summary, gene_type == "prioritized"), aes(x=stage, y=mean, group=structure_acronym, linetype="Brain regions"), color='gray') #color=structure_acronym
#p <- p + scale_colour_manual(name='BLALB', values=c('Brain regions'='grey'), guide='legend') # scale_linetype_manual
p <- p + labs(title = gene_list)
p
### Adding mean (base line)
#scale_colour_hue(h=0)
p <- p + geom_line(data=subset(df.summary.mean, gene_type == "prioritized"), aes(x=stage, y=mean1, group=1), linetype='dashed', color="#d7191c", size=1.25)
p <- p + geom_line(data=subset(df.summary.mean, gene_type == "other"), aes(x=stage, y=mean1, group=1), linetype='dashed', color='black', size=1.25)
p <- p + guides(linetype=guide_legend(keywidth = 2, keyheight = 1))
p <- p + geom_errorbar(data=subset(df.summary.mean, gene_type == "prioritized"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1, group=gene_type),color='#d7191c', width=0.2)
p <- p + geom_errorbar(data=subset(df.summary.mean, gene_type == "other"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1, group=gene_type),color='black', width=0.2)
p

### Adding GILMAN to plot
#p <- p + geom_line(data=df.gilman.summary.mean, aes(x=stage, y=mean, group = cluster, linetype=cluster, color="yellow"), size=1.25) # GILMANN cluster Ia, Ib, II
p <- p + geom_line(data=df.gilman.summary.mean.clustermerged, aes(x=stage, y=mean1, group=1, linetype="Gilman et al. cI and cII"), color="#2b83ba", size=1.25)
p <- p + geom_errorbar(data=df.gilman.summary.mean.clustermerged, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#2b83ba', width=0.2)
#p <- p + guides(linetype=guide_legend(keywidth = 2, keyheight = 1))
p


########## Adding vertical line - prenatal vs. postnatal ###########
p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
p

######### Adding x-tickmarks for stage - ##########
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
# c("s1", "s2a", "s2b", "s3a", "s3b", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11")
p + scale_x_discrete(name="", labels = stage_converter) + theme(axis.text.x = element_text(angle = 45, hjust = 1))






############################## PLAYING WITH GGPLOT LEGENDS #######################
g <- guide_legend("Legend")
p + guides(colour=g, linetype=g)

p + guides(color=guide_legend(override.aes=list(color=c("red","blue", "black"),linetype=c(1,0,1))))


p + labs(colour = "Cylinders") # The labs function also modifies legend labels


guides(fill=guide_legend(title=NULL)) # No legend

guides(colour = guide_legend(override.aes = list(size=3)))


#http://stackoverflow.com/questions/18060116/adding-legend-to-ggplot-when-lines-were-added-manually
geom_line(aes(x, y, color="My Line"), data=line.data) +
  scale_color_manual(values=c("setosa"="blue4", "versicolor"="red4",
                              "virginica"="purple4", "My Line"="gray"))






