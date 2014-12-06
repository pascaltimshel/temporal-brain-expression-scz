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

####################################### LOAD DATA ###########################################
load(file="RData/data_rnaseq_expression_processed.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt

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
################## T-test: prenatal vs postnatal ###############
### ASSUMPTION check: histogram of response variable - distribution must be "close" to normal
#ggplot(data=subset(df.expression_matrix.clean.melt, gene_type=="prioritized"), aes(x=value)) + geom_histogram(aes(y = ..density..)) + labs(title="Distribution of *selected* [priori/associated] expression values. Overlayed with dnorm(2.5,2)") + stat_function(fun = dnorm, colour = "red", arg = list(mean = 2.5, sd=2))
# other genes
##ggplot(data=subset(df.expression_matrix.clean.melt, gene_type=="other"), aes(x=value)) + geom_histogram(aes(y = ..density..)) + labs(title="Distribution of 'other' genes")


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
#save(df.fit.stage.obs, file="RData/fit_stage_obs_rnaseq_associated.RData")
#save(df.fit.stage.obs, file="RData/fit_stage_obs_rnaseq_prioritized.RData")

################## T-test: identifing most highly expressed stage - and testing for difference to other stages ###############

df.expression_matrix.clean.melt.prioritized <- subset(df.expression_matrix.clean.melt, gene_type=="prioritized")
##### WITHOUT "prior" median for structure
df.mean.prio <- ddply(df.expression_matrix.clean.melt.prioritized, c("stage"), summarise,
                               mean1 = mean(value, na.rm=TRUE),
                               sd1   = sd(value, na.rm=TRUE))
##### WITH "prior" median for structure
#df.mean.prio <- ddply(ddply(df.expression_matrix.clean.melt.prioritized, .(stage, structure_acronym), summarise, mean=median(value, na.rm=TRUE)), .(stage), summarise, mean1=mean(mean, na.rm=TRUE),  sd1=sd(mean, na.rm=TRUE))
ggplot(data=df.mean.prio, aes(x=stage, y=mean1, group=1)) + geom_line() + geom_errorbar(aes(ymin=mean1-sd1, ymax=mean1+sd1))
df.mean.prio[order(-df.mean.prio$mean),] # s11 has the highest mean expression

##### WITHOUT "prior" median for structure
# stage     mean       sd
# 12   s11 3.198934 1.711310
# 7     s6 3.164335 1.735172
# 10    s9 3.122140 1.737527

##### WITH "prior" median for structure
# stage    mean1       sd1
# 12   s11 3.344813 0.3392914
# 7     s6 3.268256 0.3266119
# 10    s9 3.172535 0.2654217

## CONCLUSION: 
# The rank of the stages does not changes with the procedure of "averaging" samples. 
# We use the "WITHOUT "prior" median for structure" approach

### TEST IT!
stage_to_test_against <- "s11"
df_test_against <- subset(df.expression_matrix.clean.melt.prioritized, stage==stage_to_test_against)
df_rest <- subset(df.expression_matrix.clean.melt.prioritized, stage!=stage_to_test_against)
ttest_stage <- function(df) {
  x <- df_test_against$value
  y <- df$value
  t.test(x, y, alternative="greater")$p.value # alternative = "greater" is the alternative that x has a larger mean than y.
  #t.test(x, y, alternative="greater") # alternative = "greater" is the alternative that x has a larger mean than y.
}
#### Making tests
#df.ttest_stage.list <- dlply(df_rest, .(stage), ttest_stage)
df.ttest_stage <- ddply(df_rest, .(stage), ttest_stage)

#### Sorting by p-value
df.ttest_stage.ordered <- df.ttest_stage[order(-df.ttest_stage[,2]),]
df.ttest_stage.ordered
#### Get stages that are NOT significant after BONFERRONI CORRECTION
df.ttest_stage.ordered[df.ttest_stage.ordered[,2] > 0.05/10,]
# stage         V1
# 7     s6 0.21841097
# 8     s7 0.02977490
# 10    s9 0.03112255
