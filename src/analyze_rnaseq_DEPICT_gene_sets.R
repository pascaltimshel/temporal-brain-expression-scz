############### SYNOPSIS ###################
# This script analyzes correlation between significant DEPICT genesets (GS_sig) and BrainSpan RNAseq data set
# Date created: 2015-02-11

### OBS ###
# The file "file.depict.gene_set" should contain ENSEBL gene IDs in the first column


### Steps ###
# 1) Match ENSEMBL gene IDs: subset DEPICT genes (~20k) to those present in BrainSpan (~50k)
# 2) Make correlation between GS_sig and
  # a) stage (12 stages; perhaps full resolutions): identify time points of active/inactive pathways
  # b) structure
# 3) Plot the correlations
  # a) stage: correlation (r) vs. time; multiple lines
  # b) structure: barplot?
# 4) Diagnositcs #1: Plot scatter plot underlying the correlation. 
  # Multipanel (facet) - one for each GS_sig
  # a) (mean) expression vs Z-score



############################################


#library(plyr)
library(dplyr)
library(tidyr)
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

str(df.expression_matrix.clean.melt)
levels(df.expression_matrix.clean.melt$stage)

########### READ prioritization file ###########
### Setting prioritization ###
file.gene_prioritization <- "../gene_lists/gene_prioritization.txt"
#file.gene_prioritization <- "../gene_lists/gene_associated.txt"
gene_list <- basename(file_path_sans_ext(file.gene_prioritization))

### Read file ###
df.gene_prioritization <- read.csv(file.gene_prioritization,h=T)
## Check mapping - compare to "df.expression_matrix.clean" NOT MOLTEN, will be too slow
sum(df.expression_matrix.clean$ensembl_gene_id %in% df.gene_prioritization[,1]) # mapped genes
sum(!df.gene_prioritization[,1] %in% df.expression_matrix.clean$ensembl_gene_id) # non-mapped genes


########### READING/LOADING DEPICT Gene Sets #########
### Read tab file - ONLY ONCE
#file.depict.gene_set <- "../data/daner_PGC_SCZ52_0513a.gz.p4.clump.areator.sorted.1mhc1_1kg_sign_genesets.tab"
#df.depict.gene_set <- read.table(file.depict.gene_set, sep="\t", h=T, stringsAsFactors=F) # first column is ENSEMBL gene IDs
#colnames(df.depict.gene_set)[1] <- "ensembl_gene_id"
#colnames(df.depict.gene_set)
#save(df.depict.gene_set, file="RData/df.depict_gene_02-17-2015.RData")

### Load RData
load(file="RData/df.depict_gene_02-17-2015.RData") # df.depict.gene_set
str(df.depict.gene_set)

### Check mapping - compare to "df.expression_matrix.clean" NOT MOLTEN, will be too slow
sum(df.expression_matrix.clean$ensembl_gene_id %in% df.depict.gene_set$ensembl_gene_id) # mapped genes
sum(!df.depict.gene_set$ensembl_gene_id %in% df.expression_matrix.clean$ensembl_gene_id) # non-mapped genes


########### Subset genes in "df.expression_matrix.clean.melt" to those in DEPICT #########
str(df.expression_matrix.clean.melt)

### Subsetting data frame
df.expression_matrix.clean.melt.sub.depict <- subset(df.expression_matrix.clean.melt, ensembl_gene_id %in% df.depict.gene_set$ensembl_gene_id)
### Converting "ensembl_gene_id" to *CHARACTER* 
df.expression_matrix.clean.melt.sub.depict$ensembl_gene_id <- as.character(df.expression_matrix.clean.melt.sub.depict$ensembl_gene_id)

### *REMOVE* variables that we do not need anymore
rm(df.expression_matrix.clean, df.expression_matrix.clean.melt)


#####################################################################################
################################ Summarizing data ###################################
#####################################################################################

### Test #1: inspecting results of grouping
#g <- group_by(df.expression_matrix.clean.melt.sub.depict, ensembl_gene_id, stage)
#g # --> same dimensions [local data frame [10,473,188 x 10] ] as "df.expression_matrix.clean.melt.sub.depict". This is likely because the data frame is already MOLTEN
   # --> you also get the same dimensions when doing "group_by(ensembl_gene_id, stage, structure_acronym)"


###### *dply* - Median per STAGE for each gene ######## 
str(df.expression_matrix.clean.melt.sub.depict)
### Run summarize
system.time(df.summary.stage <- df.expression_matrix.clean.melt.sub.depict %>%
              group_by(ensembl_gene_id, stage) %>%
              summarise(
                mean = median(value, na.rm=TRUE)
              )
)
# Time --> 19.503 (elapsed) [5 x faster! than ddply]
df.summary.stage
# Source: local data frame [239,844 x 3]

###### *dply* - Median per STRUCTURE for each gene ######## 
system.time(df.summary.structure_acronym <- df.expression_matrix.clean.melt.sub.depict %>%
              group_by(ensembl_gene_id, structure_acronym) %>%
              summarise(
                mean = median(value, na.rm=TRUE)
              )
)
df.summary.structure_acronym
# Source: local data frame [519,662 x 3]

#####################################################################################
##############################  MELT DEPICT Gene Sets ###############################
#####################################################################################

str(df.depict.gene_set)
df.depict.gene_set.melt <- melt(df.depict.gene_set, id.vars=c("ensembl_gene_id"), variable.name="gene_set", value.name="gene_set_z_score")
str(df.depict.gene_set.melt)


#####################################################################################
################################### Joining data ####################################
#####################################################################################

###### Background info - testing ########
# full_join(x, y) #--> should give two new columns with all
# inner_join(x, y) # --> same as full_join(). This is because the two data frames contain exactly the same ENSEMBL gene IDs. If df.x had more (or less) gene IDs, the resulting joined df would be smaller.
                # --> **TRY to do a inner_join with "df.expression_matrix.clean.melt" instead. Should give same result because intersection is used.
# semi_join(x, y)
# If a match is not unique, a join will add all possible combinations (the Cartesian product) of the matching observations:
# ---> use semi_join() to avoid duplications; semi_join() only ever remove observations.
# ---> BUT REMEMBER: return all rows from x where there are matching values in y, keeping just columns from x.

################ Full join ##############
### STAGE ###
#system.time(df.summary.stage.full_join <- df.summary.stage %>% full_join(df.depict.gene_set.melt, by="ensembl_gene_id"))
  # --> 234.917 s = 4 min
#save(df.summary.stage.full_join, file="RData/df.summary.stage.full_join.RData")

# Source: local data frame [34,297,692 x 5]
# Groups: ensembl_gene_id
# 
#    ensembl_gene_id stage     mean                      gene_set gene_set_z_score
# 1  ENSG00000000003   s2a 4.640608   Decreased.Vertical.Activity         1.985478
# 2  ENSG00000000003   s2a 4.640608                      Dendrite        -1.119173
# 3  ENSG00000000003   s2a 4.640608 Abnormal.Locomotor.Activation         1.444252
load(file="RData/df.summary.stage.full_join.RData") # df.summary.stage.full_join
df.summary.stage.full_join

### structure_acronym ###
#system.time(df.summary.structure_acronym.full_join <- df.summary.structure_acronym %>% full_join(df.depict.gene_set.melt, by="ensembl_gene_id"))
#df.summary.structure_acronym.full_join

#####################################################################################
######################### Correlate Genes - dply ####################################
#####################################################################################

###################### Correlation *STAGE* ###################### 
### step #1: create objects
system.time( # --> time =~ 50 s
  df.correlation.object.stage <- df.summary.stage.full_join %>%
    group_by(stage, gene_set) %>%
      ## *OBS*: cannot convert htest object to data frame: "cannot coerce class ""htest"" to a data.frame"
      ## Thus we need TWO STEPs: 1) create data.frame with htest object [this step]; 2) extract features/variables from objects [next step]
      do( # REMEMBER: NAMED arguments inside do() becomes LIST COLUMNS
        obj.htest=cor.test(.$mean, .$gene_set_z_score, na.action="na.omit"),
        # ^^ na.action="na.omit" [DEFAULT] --> will only keep PAIRWISE COMPLETE OBSERVATIONS
        
        cor = cor(.$mean, .$gene_set_z_score), # this is just for VALIDATION, i.e. check that test.cor and cor gives the same PEARSON correlation
        count = length(.$ensembl_gene_id), # cannot use: n(.)
        count_unique_genes = unique(.$ensembl_gene_id) #n_distinct(.$ensembl_gene_id)
      )
)
df.correlation.object.stage


### step #2: extract objects 
df.correlation.extract.stage <- df.correlation.object.stage %>% 
  do(data.frame(
    stage = .$stage,
    gene_set = .$gene_set,
    statistic = .$obj.htest[["statistic"]],
    p.value = .$obj.htest[["p.value"]],
    cor.estimate = .$obj.htest[["estimate"]],
    count = .$count,
    cor = .$cor
    #alternative = .$obj.htest[["alternative"]], # --> "two.sided"
    #method = .$obj.htest[["method"]], # --> "Pearson's product-moment correlation"
  ))

### Sorting
df.correlation.extract.stage <- df.correlation.extract.stage %>% arrange(desc(abs(statistic)))
df.correlation.extract.stage %>% arrange(statistic)

###################### Correlation *structure_acronym* ###################### 
### step #1: create objects 
# ...
### step #2: extract objects
# ...


#####################################################################################
######################### Line plot: correlation (r) vs time (stage) ###############
#####################################################################################
df.correlation.extract.stage
#df.correlation.extract.stage %>% slice()



### Top examples
gene_sets_of_interest <- c("Synaptosome", "TUBA1A.protein.complex", "Dendrite")
for (gs in gene_sets_of_interest) {
  p <- ggplot(df.correlation.extract.stage %>% filter(gene_set==gs), aes(x=stage, y=cor.estimate, group=1)) + geom_line() + geom_point(aes(size=log10(statistic+1)))
  #p <- p + geom_text()
  p <- p + labs(title=gs)
  filename.plot <- sprintf("tmp_figs/gene_sets_of_interest.%s.pdf", gs)
  ggsave(file=filename.plot)
}

#### Getting Gene Sets with correlations at any given time point MORE extreme (high/low) than a given threshold ####
### Active (r>=0.3)
gs_active <- df.correlation.extract.stage %>%
  group_by(gene_set) %>%
  summarise(
    cor.max = max(cor),
    cor.mean = mean(cor)
  ) %>%
  filter(cor.max >= 0.3) %>%
  arrange(desc(cor.max))
gs_active

### Inactive (r<=-0.2)
gs_inactive <- df.correlation.extract.stage %>% 
  group_by(gene_set) %>%
  summarise(
    cor.min = min(cor),
    cor.mean = mean(cor)
  ) %>%
  filter(cor.min <= -0.2) %>%
  arrange(cor.min)
gs_inactive


### Overlap
intersect(gs_active$gene_set, gs_inactive$gene_set) # --> nothing
cat(as.character(gs_active$gene_set), sep="\n")

######################### Plot most active/inactive pathways #########################

#x %>% f(y) -> f(x,y)

df.correlation.extract.stage %>% 
  mutate(
    plot_previously = gene_set %in% gs_active$gene_set
  ) %>%
  group_by(gene_set) %>%
  summarise(
      increasing = ifelse(stage=="")
  )

p <- ggplot(df.correlation.extract.stage, aes(x=stage, y=cor.estimate)) + geom_boxplot()
p + labs(title="SCZ")
ggsave(file="tmp_figs/boxplot_scz.pdf", w=8, h=6)

### New - color by plot_previously
p <- ggplot(df.correlation.extract.stage %>% mutate(plot_previously = gene_set %in% gs_active$gene_set), aes(x=stage, y=cor.estimate, group=gene_set)) + geom_line(aes(color=plot_previously)) + geom_point(aes(size=log10(abs(statistic)+1)))
p
## new - all
#p <- ggplot(df.correlation.extract.stage, aes(x=stage, y=cor.estimate, group=gene_set)) + geom_line(aes(color=gene_set)) + geom_point(aes(size=log10(abs(statistic)+1)))
p + guides(color=FALSE) + geom_hline(y=0.3)
ggsave("tmp_figs/tmp1.pdf", w=20, h=20)

### Active
p <- ggplot(df.correlation.extract.stage %>% filter(gene_set %in% gs_active$gene_set), aes(x=stage, y=cor.estimate, group=gene_set)) + geom_line(aes(color=gene_set)) + geom_point(aes(size=log10(statistic+1)))
p + guides(color=FALSE)
#p + theme(legend.position="bottom") + guides(color=guide_legend(nrow=6, ncol=4))
#ggsave(file="tmp_figs/active_pathways_n24_no-legend.pdf")
#ggsave(file="tmp_figs/active_pathways_24.pdf")

### Inactive
x <- df.correlation.extract.stage %>% filter(gene_set %in% gs_inactive$gene_set)
p <- ggplot(df.correlation.extract.stage %>% filter(gene_set %in% gs_inactive$gene_set), aes(x=stage, y=cor.estimate, group=gene_set)) + geom_line(aes(color=gene_set)) + geom_point(aes(size=log10(abs(statistic)+1)))
p
#ggsave(file="tmp_figs/inactive_pathways_24.pdf")


#####################################################################################
######################### SCATTER PLOT: expr vs z_score #################################
#####################################################################################

### Get a few gene sets to make scatter plot
n_gene_sets <- 2
#gs2plot <- sample_n(gs_active, n_gene_sets)$gene_set # random sampling
gs2plot.active <- slice(gs_active, 1:n_gene_sets)$gene_set # pick top ranking [REMEMBER the data frame is sorted!]
gs2plot.inactive <- slice(gs_inactive, 1:n_gene_sets)$gene_set
gs2plot.active
gs2plot.inactive

### Active
df.summary.stage.full_join %>% filter(gene_set %in% gs2plot.active)
p <- ggplot(df.summary.stage.full_join %>% filter(gene_set %in% gs2plot.active), aes(x=mean, y=gene_set_z_score)) + geom_point()
p + facet_grid(gene_set ~ stage) # facet_grid(y ~ x) 
filename.plot <- sprintf("tmp_figs/facet_grid_active_top%s.pdf", n_gene_sets)
ggsave(file=filename.plot, width=20, height=20)

### Active - WITH COLORING OF PRIORITIZED GENES
df.x <- df.summary.stage.full_join %>% 
  filter(gene_set %in% gs2plot.active) %>%
  mutate(
    prioritized = ensembl_gene_id %in% df.gene_prioritization[,1]
  )
p <- ggplot(df.x, aes(x=mean, y=gene_set_z_score)) + geom_point(aes(color=prioritized))
p + facet_grid(gene_set ~ stage) # facet_grid(y ~ x) 
filename.plot <- sprintf("tmp_figs/facet_grid_active_top%s_color_prioritized.pdf", n_gene_sets)
ggsave(file=filename.plot, width=20, height=20)

#### TODO:
# 1) add correlation TEXT to each subplot/grid.





#####################################################################################
#################################### GARBAGE ########################################
#####################################################################################

gs_active <- df.correlation.extract.stage %>% 
  filter(cor.estimate >= 0.3) %>%
  select(gene_set) %>%
  distinct() 
#%>% 
#unlist()
#as.data.frame()
#as.character()
gs_active$gene_set