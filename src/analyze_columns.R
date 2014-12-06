############### SYNOPSIS ###################
# This script will analyze the columns file
# The script will produce the following
# 1) CONTINGENCY TABLE plot
# 2) Marginal barplot
# The script also performs comparison of the samples for Microarray and RNAseq
############################################


library(ggplot2)
library(reshape2)
library(plyr)

rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

########### SOURCING ###########
source("function_def_stages.R", echo=TRUE)

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


###############################


#### MICROARRAY
#file.columns <- "../data/141031/microarray/columns_metadata.csv"
#### RNAseq
file.columns <- "../data/141031/rnaseq/columns_metadata.csv"

########### READ columns file ###########
df.columns <- read.csv(file.columns,h=T,row.names=1)
#### Add stage column
df.columns$stage <- as.factor(sapply(df.columns$age, function(x) {names(stages)[sapply(stages, function(stage) x %in% stage)]}))
#### Sort factor levels of "stage"
df.columns$stage <- factor(df.columns$stage, levels(df.columns$stage)[match(order.stages, levels(df.columns$stage))])


######### Number of factor levels for variables
paste("Number of unique donors:", length(unique(df.columns$donor_id)))
paste("Number of unique age:", length(unique(df.columns$age)))
paste("Number of unique structure_acronum:", length(unique(df.columns$structure_acronym)))


############ Plot contingency table
table.contingency <- table(df.columns$stage, df.columns$structure_acronym)
df.contingency.melt <- melt(table.contingency, varnames=c("stage", "structure"), value.name=c("count"))
p <- ggplot(df.contingency.melt, aes(x=stage,y=structure)) + geom_tile(aes(fill=count))
p <- p + geom_text(aes(label=count), size=rel(5), color="white", family="mono")
p

### Margins
df.margin.stage <- data.frame(margin.table(table.contingency, 1))
df.margin.structure <- data.frame(margin.table(table.contingency, 2))

###### Stage plot
p <- ggplot(df.margin.stage, aes(x=Var1, y=Freq)) + geom_bar(stat='identity') + labs(x="Stage", y="#Samples")
p
p <- do_stage_converter(p)
p + theme(text = element_text(size=15))

##### "Age" plot | using df.columns
df.columns.age <- as.data.frame(table(df.columns$age))
variable_split <- strsplit(as.character(df.columns.age$Var1), " ") # list of list with two elements
df.columns.age$value <- as.numeric(sapply(variable_split, "[[", 1)) # integer/numeric
df.columns.age$unit <- sapply(variable_split, "[[", 2) # unit {pcw, mos, yrs}
str(df.columns.age)
# sorting data frame
df.columns.age <- df.columns.age[with(df.columns.age, order(match(unit, c("pcw", "mos", "yrs")), value)), ]
str(df.columns.age)
df.columns.age$Var1
# converting "age"/Var1 to factor and sorting
df.columns.age$Var1 <- factor(df.columns.age$Var1, levels=df.columns.age$Var1) # OBS: ordered=T
levels(df.columns.age$Var1)
str(df.columns.age)
## PLOT IT
p <- ggplot(df.columns.age, aes(x=Var1, y=Freq)) + geom_bar(stat='identity') + labs(x="Age", y="#Samples")
p + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=rel(1.15)))

##### Structure plot
ggplot(df.margin.structure, aes(x=Var1, y=Freq)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 35, hjust = 1)) + theme(text = element_text(size=15)) + labs(x="Structure", y="#Samples")





# COMPARISON OF MICROARRAY VS RNASEQ  ----------------------------------------------------------------------

#### MICROARRAY
file.columns.marray <- "/Users/pascaltimshel/p_scz/brainspan/data/141031/microarray/columns_metadata.csv"
#### RNAseq
file.columns.rnaseq <- "/Users/pascaltimshel/p_scz/brainspan/data/141031/rnaseq/columns_metadata.csv"

########### READ columns file ###########
df.columns.marray <- read.csv(file.columns.marray,h=T,row.names=1)
df.columns.rnaseq <- read.csv(file.columns.rnaseq,h=T,row.names=1)

########## Finding overlapping observations

df <- ldply(.data=list(marray=df.columns.marray, rnaseq=df.columns.rnaseq))
bool.dup <- duplicated(df[,-1]) | duplicated(df[,-1], fromLast = TRUE) # do not include ".id" column when considering duplicates
df.dup <- df[bool.dup,]
df.no.dup <- df[!bool.dup,] # observations/columns that are not duplicated!

df.no.dup.stats <- ddply(df.no.dup, .(.id), summarise, 
                         n=length(.id))

#marray=59 (unique samples not in rnaseq)
#rnaseq=91 (unique samples not in marray)
#diff=91-59=32

######### Number of donors
length(unique(df.columns.marray$donor_id)) #---> 35 donors
length(unique(df.columns.rnaseq$donor_id)) #---> 42 donors

########## Donor overlap
length(intersect(df.columns.marray$donor_id, df.columns.rnaseq$donor_id)) # 35. There are 35 unique donor_id in marray
all(df.columns.marray$donor_id %in% df.columns.rnaseq$donor_id) # TRUE --> marray donors is a subset of rnaseq
df.news.rnaseq.donor <- df.columns.rnaseq[!(df.columns.rnaseq$donor_id %in% df.columns.marray$donor_id), ]
df.news.rnaseq.donor
### New donors
cat(unique(as.character(df.news.rnaseq.donor$donor_id)), sep='; ') # 12295; 263195015; 12889; 12890; 12981; 12289; 12832
### Age of new samples
table(as.character(df.news.rnaseq.donor$age))
cat(unique(as.character(df.news.rnaseq.donor$age)), sep='; ') #35 pcw; 37 pcw; 4 mos; 8 yrs; 11 yrs; 19 yrs --> "4 mos; 8 yrs;" are seen in marray

######### Age
length(intersect(df.columns.marray$age, df.columns.rnaseq$age)) # 27. marray unique age: 27
all(df.columns.marray$age %in% df.columns.rnaseq$age) # TRUE
# New age in rnaseq
cat(unique(as.character(df.columns.rnaseq[!(df.columns.rnaseq$age %in% df.columns.marray$age), "age"])), sep='; ')
# ---> 35 pcw; 37 pcw; 11 yrs; 19 yrs
?intersect
df.columns.marray$donor_id, df.columns.rnaseq$donor_id

df.columns.marray$donor_id, df.columns.rnaseq$donor_id




