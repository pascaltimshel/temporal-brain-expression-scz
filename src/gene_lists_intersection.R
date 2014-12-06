############### SYNAPSIS ###################
# Name: graphics_main1_median.R
# This script is a newer version of "graphics_main1.R".
# The script LOADs EITHER microarray data OR PROCESSED RNAseq data.
# It was created to be able to handle multiple gene list for plotting. The script uses plyr to load in multiple gene lists.
# You may choose what files to load
# The script also process and loads GILMAN gene clusters (1, 1a, 1b and 2)
# The script was used to generate PUBLICATION READY GRAPHICS
# THIS VERSION TAKES THE MEDIAN GENE EXPRESSION VALUES.
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/p_scz/brainspan/src"
setwd(wd)

############################# READING GENE LISTs #################################
path.datafiles <- '/Users/pascaltimshel/p_scz/brainspan/gene_lists'

###### Read into a list of files - PATTERN VERSION - read ALL .txt files in directory:
#files <- list.files(path = path.datafiles, pattern = "*.txt", full.names = TRUE) #full path
#names(files) <- list.files(path = path.datafiles, pattern = "*.txt") # filename
#cat(names(files), sep="\n")

###### Read SPECIFIC FILES:
filenames2read <- c("gene_prioritization.txt", "gene_associated.txt", "gene_nearest.txt")
#filenames2read <- c(filenames2read, "gene_psd_human.txt", "gene_psd_mouse.txt")
filenames2read <- c(filenames2read, "gilman_nn_2012_cluster1.ens", "gilman_nn_2012_cluster1a.ens", "gilman_nn_2012_cluster1b.ens", "gilman_nn_2012_cluster2.ens")
filenames2read <- c(filenames2read, "gulsuner_S3A_damaging_cases.ens")
files <- as.list(paste(path.datafiles, filenames2read, sep="/"))
names(files) <- filenames2read
files
list_of_data <- llply(files, read.csv)#row.names = 1 --> NO!, stringsAsFactors = FALSE

names(list_of_data)



############## MAKING HEATMAP PLOT ###################
dflist <- list_of_data # copying dflist.new [TODO: make the below code a function]
mat.intersect <- matrix(data=NA, nrow=length(dflist), ncol=length(dflist), dimnames=list(names(dflist), names(dflist)))
for (i in 1:length(dflist)) {
  #gene_list = dflist[[i]][,1]
  for (j in i:length(dflist)) {
    #cat("i=",i,"j=",j,"\n")
    mat.intersect[i,j] <- length(intersect(dflist[[i]][,1], dflist[[j]][,1])) # Symmetric matrix
    mat.intersect[j,i] <- length(intersect(dflist[[i]][,1], dflist[[j]][,1])) # Symmetric matrix
    # TRY ODDs ratio
    #mat.intersect[i,j] <- length(intersect(dflist[[i]][,1], dflist[[j]][,1]))/length(dflist[[j]][,1]) # NORMALIZING by "column" length
  }
}

melt.mat.intersect <- melt(mat.intersect)

### Renaming gene lists
rename_map <- c("gene_associated.txt"="Associated Genes", 
                "gene_nearest.txt"="Nearest Genes", 
                "gene_prioritization.txt"="Prioritized Genes",
                "gilman_nn_2012_cluster1.ens"="Gilman et al. cluster I",
                "gilman_nn_2012_cluster1a.ens"="Gilman et al. cluster Ia",
                "gilman_nn_2012_cluster1b.ens"="Gilman et al. cluster Ib",
                "gilman_nn_2012_cluster2.ens"="Gilman et al. cluster II",
                "gulsuner_S3A_damaging_cases.ens"="Gulsuner et al. table S3A",
                "gene_psd_human.txt"="Post Synaptic Genes (Human)", 
                "gene_psd_mouse.txt"="Post Synaptic Genes (Mouse)")
levels(melt.mat.intersect$Var1)
melt.mat.intersect$Var1 <- revalue(melt.mat.intersect$Var1, rename_map)
melt.mat.intersect$Var2 <- revalue(melt.mat.intersect$Var2, rename_map)
levels(melt.mat.intersect$Var1)


q1 <- ggplot(melt.mat.intersect, aes(Var1, Var2, fill = value)) + geom_tile()
base_size <- 12
q1 <- q1 + theme_grey(base_size = base_size) + labs(x = "", y = "", title="Gene list intersections") 
q1 <- q1 + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
#q1 <- q1 + theme(axis.text.x=element_text(size = base_size*0.8, angle=330, hjust=0, colour="grey50"))
q1 <- q1 + theme(axis.text.x=element_text(size = base_size*0.8, angle=315, hjust=0, colour="grey50"))
q1 <- q1 + geom_text(aes(label = value), size=4)
#q1 <- q1 + coord_fixed() # equal axis
q1 + scale_fill_gradient(low = "#fee8c8",  high = "#e34a33") # light_orange-->orange
#q1 + scale_fill_gradient(low = "#ece2f0",  high = "#1c9099") # baise-->green
#q1 + scale_fill_gradient(low = "#f0f0f0",  high = "#636363") # black-->white
#q1 + scale_fill_gradient(low = "blue",  high = "green")

#gene_list_intersection_c1-8x6.pdf


#q1 + scale_color_brewer()
# http://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/


#################################### Make VennDiagram ###############################



