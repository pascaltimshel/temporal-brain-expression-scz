############### SYNOPSIS ###################
# This is a helper function for another script
############################################


library(ggplot2)
library(reshape2)
library(plyr)



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
# p <- ggplot(df.columns.age, aes(x=Var1, y=Freq/max(Freq))) + geom_bar(stat='identity') + labs(x="Age", y="#Samples")
# p + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=rel(1.15)))
