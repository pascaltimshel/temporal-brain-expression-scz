library(plyr)
library(ggplot2)
library(reshape2)

rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

################################# LOAD *UNPROCESSED* DATA ###################################

#load(file="RData/data_rnaseq_expression_unprocessed_full_res.RData")
#df.expression_matrix.clean <- df.expression_matrix.clean.unprocessed
#df.expression_matrix.clean.melt <- df.expression_matrix.clean.melt.unprocessed

################################# PROCESSING DATA #####################################
########## Quantiles and values > W_param RPKM
W <- 200
sum(df.expression_matrix.clean.melt$value > W)
# W(50)=542737
# W(100)=233217
genes_affected_by_winsorizing <- unique(df.expression_matrix.clean.melt[df.expression_matrix.clean.melt$value > W, c("ensembl_gene_id")])
length(genes_affected_by_winsorizing) 
# W(50)=4586
# W(100)=2349

######### Plot quantiles
step_size <- 0.005 # the step size defines the "fine grainedness" of the plot.
t <- quantile(df.expression_matrix.clean.melt$value, probs=seq(0,1,step_size))
t # named vector
###  CONTRUCT data frame for lines: get the values of the quantiles of the interest
marks_probs <- c(0.95, 0.98, 0.99, 0.995)
marks_value <- quantile(df.expression_matrix.clean.melt$value, probs=marks_probs)
df.marks <- data.frame(marks_probs, marks_value, marks_text=names(marks_value))
### Making plot
quantile_step_vec <- seq(0, 1, step_size) # this will be the y-axis
p <- ggplot(data=data.frame(quantile_step_vec, t), aes(x=t, y=quantile_step_vec))
#p <- p + geom_step() # this is nice too
p <- p + geom_point()
p <- p + xlim(0,200) # removes data
# HORIZONTAL LINE
p <- p + geom_segment(data=df.marks, aes(x=0, xend=marks_value, y=marks_probs, yend=marks_probs), linetype="dashed") # x=-Inf would also work
# VERTICAL LINE
p <- p + geom_segment(data=df.marks, aes(x=marks_value, xend=marks_value, y=0, yend=marks_probs), linetype="dashed") # y=-Inf would also work
# ANNOTATE POINTS
p <- p + geom_text(data=df.marks, aes(x=marks_value, y=marks_probs, label=marks_text), vjust=1.1, hjust=-0.1) # (x,y) points
p <- p + geom_text(data=df.marks, aes(x=marks_value, y=0, label=round(df.marks$marks_value, 1)), vjust=-0.1, hjust=-0.1) # (x, y=0) points
p + labs(x="RPKM", y="quantile", title="RNAseq winsorizing effect visualization")
#geom_segment(aes(x=0, y=0.99, xend=t, yend=0.99))


########## ALTERNATIVES 
### Option 1
#Empirical Cumulative Density Function
#ggplot(df.expression_matrix.clean.melt, aes(x=value)) + stat_ecdf(n=5) # I think this will just hang for a long time

### Option 2
#ggplot(df, aes(x=1:5, y=cumsum(val))) + geom_line() + geom_point()

### Option 3
# geom_step(mapping = NULL, data = NULL, stat = "identity", position = "identity",)
# Connect observations by stairs.
