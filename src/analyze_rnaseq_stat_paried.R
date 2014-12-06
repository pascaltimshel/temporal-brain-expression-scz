############### SYNOPSIS ###################
# This script does statistical testing on the RNAseq data set
# This script is a more SIMPEL (stripped) version of "analyze_marray_stat.R"
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

####################################### SOURCE ###########################################
source("multiplot.R")

####################################### LOAD DATA ###########################################
load(file="RData/data_rnaseq_expression_processed_full_res.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt

########### Setting prioritization #########
file.gene_prioritization <- "../gene_lists/gene_prioritization.txt"
#file.gene_prioritization <- "../gene_lists/gene_associated.txt"

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


################################# STATISTICAL TESTs ##################################
################## PREPARING FOR TEST ###############
## subsetting on prioritized/associated genes
df.expression_matrix.clean.melt.priori <- subset(df.expression_matrix.clean.melt, gene_type=="prioritized")

## calculating mean for each gene in prenatal and postnatal
df.natal.gene.mean <- ddply(df.expression_matrix.clean.melt.priori, .(ensembl_gene_id, natal), summarise,
                            mean=mean(value, na.rm=T))


################## T-test: prenatal vs postnatal | *PAIRED* ###############
# df.tmp.prenatal <- subset(df.natal.gene.mean, natal=="prenatal")
# df.tmp.postnatal <- subset(df.natal.gene.mean, natal=="postnatal")
# fit.test <- t.test(df.tmp.postnatal$mean, df.tmp.prenatal$mean, alternative="greater", paired=TRUE)
# fit.test # USING VECTORS ---> 0.3211001

### t.test
levels(df.natal.gene.mean$natal)
fit1 <- t.test(mean~natal, data=df.natal.gene.mean, alternative="greater", paired=TRUE)
fit1
fit1$p.value
### calculation odd's ratio
df.or <- ddply(df.natal.gene.mean, .(natal), summarise,
               group_mean=mean(mean, na.rm=T))
df.or
or_log2_corrected <- (2^df.or[df.or$natal=="postnatal","group_mean"]-1)/(2^df.or[df.or$natal=="prenatal","group_mean"]-1)
#formula: or = ((2^mu1)-1)/((2^mu2)-1)
or_log2_corrected

### plotting distribution of mean values
q1 <- ggplot(df.natal.gene.mean, aes(x=mean, fill=natal)) + geom_histogram(aes(y=..density..)) +geom_density(alpha=0.8)
### plotting box plot
q2 <- ggplot(df.natal.gene.mean, aes(x=natal, y=mean, fill=natal)) + geom_boxplot()
multiplot(q1,q2)

#######################################################################################################################
################## T-test: identifing most highly expressed stage - and testing for difference to other stages ###############
################# PREPARING FOR TEST - step#1: calculate mean of each gene for each stage ###############
##### WITHOUT "prior" median for structure  ######
# ---> We do not use approach because it is NOT similar to how we plot data
# ---> This approach makes sense since it gives a "weighted" average 
#      ----> here structures with low samples does NOT count as much as structures with many samples
#df.mean.prio <- ddply(df.expression_matrix.clean.melt.priori, .(ensembl_gene_id, stage), summarise,
#                            mean1=mean(value, na.rm=T))

##### WITH "prior" median for structure ######
# ---> We use this approach because it IS SIMILAR TO HOW WE PLOT DATA 
# ---> This is not an weighted average approach
df.mean.prio <- ddply(ddply(df.expression_matrix.clean.melt.priori, .(ensembl_gene_id, stage, structure_acronym), summarise, mean=median(value, na.rm=TRUE)), .(ensembl_gene_id, stage), summarise, mean1=mean(mean, na.rm=TRUE),  sd1=sd(mean, na.rm=TRUE))
################# PREPARING FOR TEST - step#2: calculate mean value across all genes for each stage  ###############
### This data frame is only used for identifying the HIGHEST expressed stages and PLOTTING
df.mean.prio.stage <- ddply(df.mean.prio, .(stage), summarise,
                            mean_stage=mean(mean1),
                            sd_stage=sd(mean1))
################# PREPARING FOR TEST - step#3: plot it! #################
p <- ggplot(data=df.mean.prio, aes(x=stage, y=mean1)) + geom_line(aes(group=ensembl_gene_id, color=ensembl_gene_id))
p <- p + geom_line(data=df.mean.prio.stage, aes(x=stage, y=mean_stage, group=1), size=1.5, color="black")
p <- p + geom_errorbar(data=df.mean.prio.stage, aes(x=stage, y=mean_stage, ymin=mean_stage-sd_stage, ymax=mean_stage+sd_stage, group=1), size=1.5, color="black")
p
################# PREPARING FOR TEST - step#4: identify highest expressed stage #################
df.mean.prio.stage[order(-df.mean.prio.stage$mean_stage),] 
# ---> s11 has the highest mean expression
# stage mean_stage sd_stage
# 12   s11   3.214713 1.618984
# 7     s6   3.196868 1.638089
# 11   s10   3.177876 1.640300


################# MAKING TEST - step#5 #################
stage_to_test_against <- "s11"
df_test_against <- subset(df.mean.prio, stage==stage_to_test_against)
df_rest <- subset(df.mean.prio, stage!=stage_to_test_against)
ttest_stage <- function(df) {
  x <- df_test_against$mean1 # *OBS*: mean1
  y <- df$mean1 # *OBS*: mean1
  ### *PAIRED t-test* ###
  t.test(x, y, alternative="greater", paired=TRUE)$p.value # alternative = "greater" is the alternative that x has a larger mean than y.
  #t.test(x, y, alternative="greater", paired=TRUE) # ---> use this in COMBINATION WITH dlply
}
#### Making tests
#df.ttest_stage.list <- dlply(df_rest, .(stage), ttest_stage)
df.ttest_stage <- ddply(df_rest, .(stage), ttest_stage)

#### Sorting by p-value
df.ttest_stage.ordered <- df.ttest_stage[order(-df.ttest_stage[,2]),]
df.ttest_stage.ordered
#### Get stages that are NOT significant after BONFERRONI CORRECTION
df.ttest_stage.ordered[df.ttest_stage.ordered[,2] > 0.05/10,]
# stage          V1
# 7     s6 0.389784242
# 11   s10 0.070945047
# 5     s4 0.036666238
# 8     s7 0.012373153
# 4    s3b 0.008734785
# 10    s9 0.007094453


