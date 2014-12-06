############### SYNOPSIS ###################
# Perform cluster analysis
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

####################################### LOAD DATA ###########################################
load(file="RData/data_rnaseq_expression_processed_full_res.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt

####################################### SOURCE ###########################################
source("multiplot.R")
######## Defining stages #########
source("function_def_stages.R", echo=TRUE)
source("function_load_gene_list_data_individual_gene_plots.R", echo=TRUE)
# --> do_stage_converter
# --> df.summary
# --> df.summary.sem
# --> df.all.sem


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

### Setting priorizied factor for UN-MOLTEN df + MORE!
rownames(df.expression_matrix.clean) <- as.character(df.expression_matrix.clean$ensembl_gene_id)
df.expression_matrix.clean$gene_type <- as.factor(ifelse(df.expression_matrix.clean$ensembl_gene_id %in% df.gene_prioritization[,1], "prioritized", "other"))
table(df.expression_matrix.clean$gene_type)

### Setting priorizied factor for MOLTEN df
#df.expression_matrix.clean.melt$gene_type <- as.factor(ifelse(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.gene_prioritization[,1], "prioritized", "other"))
#table(df.expression_matrix.clean.melt$gene_type)
#str(df.expression_matrix.clean.melt)

################## subsetting on prioritized/associated genes ###############
## un-molten data frame
str(df.expression_matrix.clean, list.len=Inf)
df.expression_matrix.clean.priori <- subset(df.expression_matrix.clean, gene_type=="prioritized", select=c(-ensembl_gene_id, -gene_length, -gene_type))
## un-molten data frame
#df.expression_matrix.clean.melt.priori <- subset(df.expression_matrix.clean.melt, gene_type=="prioritized")


#http://www.statmethods.net/advstats/cluster.html
###################################### AGLOMETRIC ###########################################
# Ward Hierarchical Clustering
mydata <- df.expression_matrix.clean.priori
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

###################################### KMEANS ###########################################
# Determine number of clusters
mydata <- df.expression_matrix.clean.priori
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 3) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
kmean_res <- data.frame(ensembl_gene_id=rownames(mydata), cluster=fit$cluster)


###############################################################################

########## Read ENSEMBL-2-HGNC mapping file for prioritized genes ##############
file.map.prioritized <- '/Users/pascaltimshel/p_scz/brainspan/gene_lists/gene_prioritization.csv'
df.map.prioritized <- read.csv(file.map.prioritized, h=F)




### Subsetting genes
gene_ens <- as.character(df.map.prioritized[,1]) # ALL GENES!
df.sub <- subset(df.expression_matrix.clean.melt, (ensembl_gene_id %in% gene_ens & structure_acronym %in% c("DFC", "MFC", "OFC", "VFC")))
### constructing label: mapping ensembl-->hgnc
row_idx <- match(df.sub$ensembl_gene_id, df.map.prioritized[,1]) # No need for as.character(): match() --> factors are converted to character vectors,
hgnc_symbol <- as.character(df.map.prioritized[row_idx,2])
# replacing unmapped genes (hgnc=="") with their ENSEMBL names:
df.sub$facet_label <- ifelse(hgnc_symbol!="", hgnc_symbol, as.character(df.sub$ensembl_gene_id))

### ASSIGNING CLUSTERS TO GENES ###
match_idx <- match(df.sub$ensembl_gene_id, kmean_res$ensembl_gene_id)
df.sub$cluster <- kmean_res$cluster[match_idx]




### Making plot
save_images <- FALSE
for (i in unique(df.sub$cluster)) {
  df.sub <- df.sub[df.sub$cluster==i, ]
  n_cluster_members <- length(unique(df.sub$ensembl_gene_id))
  print(n_cluster_members)
  
  if (save_images) {
    p <- plot_cluster_wrap(df.sub)
    file_name <- paste("cluster_", i, ".pdf", sep="") # CONSIDER USING .SVG FORMAT
    pdf(file_name, width=10, height=6)
    print(p)
    dev.off()
  }
}

plot_cluster_wrap <- function(df.sub) {
  p <- ggplot(data=df.sub)
  p <- p + geom_point(data=df.sub, aes(x=stage, y=value, color=structure_acronym))
  p <- p + geom_smooth(data=df.sub, aes(x=stage, y=value, group=structure_acronym, color=structure_acronym), method="loess", alpha=0.1)
  ## add mean of prioritized
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Mean all structures"), linetype='solid', size=1)
  ## add mean of all genes
  #p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1), linetype='solid', color="Black", size=1)
  #### SETTING LEGEND
  cmap <- c("Mean all structures"="gray",
            "DFC"="#F8766D",
            "MFC"="#7CAE00",
            "OFC"="#00BFC4",
            "VFC"="#C77CFF",
            guide='legend')
  p <- p + scale_color_manual(name="Brain structure", values=cmap)
  p <- do_stage_converter(p)
  ## Setting x-axis text size (DIFFICULT!)
  p <- p + theme(axis.text.x=element_text(size=8))
  p <- p + labs(y="Mean brain expression")
  ### Setting text size for label/strip.text in facets
  p <- p + theme(strip.text.x=element_text(size=16, colour="black", face="italic"))
  ### Setting axis.title size for y-axis
  p <- p + theme(axis.title.y=element_text(size=20))
  p <- p + facet_wrap(~ facet_label)
  return(p)
}
#prioritized_genes_structure=FC_with_mean-20x20.pdf


















