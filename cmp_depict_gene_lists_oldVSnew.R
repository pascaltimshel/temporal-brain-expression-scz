library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)


#################################### ASSOCIATED ##########################################
file.a.1 <- "../gene_lists/gene_associated.txt"
file.a.2 <- "../gene_lists/depict_broad_analysis_OUTDATED_per_11-13-2014/gene_associated.txt"

df.a.1 <- read.csv(file.a.1)
df.a.2 <- read.csv(file.a.2)
x<-intersect(df.a.1[,1],df.a.2[,1])
length(x)

#################################### PRIORITIZED #########################################

file.p.1 <- "../gene_lists/gene_prioritization.txt"
file.p.2 <- "../gene_lists/depict_broad_analysis_OUTDATED_per_11-13-2014/gene_prioritization.txt"

df.p.1 <- read.csv(file.p.1)
df.p.2 <- read.csv(file.p.2)
x<-intersect(df.p.1[,1],df.p.2[,1])
length(x)

#################################### NEAREST #########################################


file.p.1 <- "../gene_lists/gene_nearest.txt"
file.p.2 <- "../gene_lists/depict_broad_analysis_OUTDATED_per_11-13-2014/gene_nearest.txt"

df.p.1 <- read.csv(file.p.1)
df.p.2 <- read.csv(file.p.2)
x<-intersect(df.p.1[,1],df.p.2[,1])
length(x)
