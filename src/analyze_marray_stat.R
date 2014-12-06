############### SYNOPSIS ###################
# This script does statistical testing on the microarray data set
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


########### Setting prioritization #########
file.gene_prioritization <- "../gene_lists/gene_prioritization.txt"
#file.gene_prioritization <- "../gene_lists/gene_associated.txt"
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


################################# STATISTICAL TESTs ##################################
################## T-test: prenatal vs postnatal ###############
### ASSUMPTION check: histogram of response variable - distribution must be "close" to normal
ggplot(data=subset(df.expression_matrix.clean.melt, gene_type=="prioritized"), aes(x=value)) + geom_histogram(aes(y = ..density..)) + labs(title="Distribution of prioritized genes expression values. Overlayed with dnorm(9,2)") + stat_function(fun = dnorm, colour = "red", arg = list(mean = 9, sd=2))

### t.test
fit1 <- t.test(value~natal, data=df.expression_matrix.clean.melt, alternative="greater", subset=df.expression_matrix.clean.melt$gene_type=="prioritized")
fit1
fit1$p.value
fit1$estimate[1]/fit1$estimate[2]
fit1$estimate[1]-fit1$estimate[2]

################## T-test: higher expression ###############
### BOXPLOT # ** takes ~10 sec to plot
ggplot(df.expression_matrix.clean.melt, aes(x=gene_type, y=value)) + geom_boxplot(aes(fill=gene_type))

### t.test ** runtime ~10s
fit.expr.higher <- t.test(value~gene_type, data=df.expression_matrix.clean.melt, alternative="less") # other < prioritized ("o" comes before "p")
fit.expr.higher
fit.expr.higher$p.value


################## T-test: test for each specific time point ###############

linmod <- function(df) {
  t.test(value~gene_type, data=df, alternative="less") # other < prioritized ("o" comes before "p")
}
list.fit.stage <- dlply(df.expression_matrix.clean.melt, .(stage), linmod) # empty factor levels are dropped by default.
df.fit.stage.obs<-ldply(list.fit.stage, function(x) {data.frame(t_statistic=x$statistic)})
qplot(stage,-t_statistic, geom="line", group=1, data=df.fit.stage.obs) + labs(title="Micro-array", y="-t statistic")

######## Save data ########
save(df.fit.stage.obs, file="RData/fit_stage_obs_marray_associated.RData")
#save(df.fit.stage.obs, file="RData/fit_stage_obs_marray_prioritized.RData")




################################# STATISTIC MODELS #######################################
head(df.expression_matrix.clean.melt)
fit <- lm(value ~ gene_type+structure_acronym+gene_type:structure_acronym, data=df.expression_matrix.clean.melt)
summary(fit)


################################# GARBAGE #######################################

