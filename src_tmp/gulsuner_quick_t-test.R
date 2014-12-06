############### SYNOPSIS ###################
## TMP TMP TMP
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

####################################### LOAD DATA ###########################################
load(file="RData/data_rnaseq_expression_processed.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt

########### Setting prioritization #########
file.gene_prioritization <- "../gene_lists/gulsuner_S3A_damaging_cases.ens"

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
################## T-test: prenatal vs postnatal - ALL STRUCTURES ###############
### t.test | postnatal > prenatal
fit1 <- t.test(value~natal, data=df.expression_matrix.clean.melt, alternative="greater", subset=df.expression_matrix.clean.melt$gene_type=="prioritized")
fit1
fit1$p.value
fit1$estimate[1]/fit1$estimate[2]
fit1$estimate[1]-fit1$estimate[2]

### t.test | postnatal < prenatal
fit1 <- t.test(value~natal, data=df.expression_matrix.clean.melt, alternative="less", subset=df.expression_matrix.clean.melt$gene_type=="prioritized")
fit1
fit1$p.value
fit1$estimate[1]/fit1$estimate[2]
fit1$estimate[1]-fit1$estimate[2]

################## T-test: prenatal vs postnatal - DFC STRUCTURE ###############
str(df.expression_matrix.clean.melt)
### t.test | postnatal < prenatal
df.sub <- subset(df.expression_matrix.clean.melt, (gene_type=="prioritized" & structure_acronym=="DFC"))
fit1 <- t.test(value~natal, data=df.sub, alternative="less")
fit1
fit1$p.value
fit1$estimate[1]/fit1$estimate[2]
fit1$estimate[1]-fit1$estimate[2]

