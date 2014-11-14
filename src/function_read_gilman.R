############### SYNAPSIS ###################
# This function is descrepated but works. Instead, use the newer "graphics" scripts to load in the Gilman clusters via plyr
# OUTPUT: df.gilman.summary.mean
############################################


####################### GILMAN LIST ################################
df.gilman.list <- list()
### 
#df.gilman.list[["cluster1"]] <- read.csv("../gene_lists/gilman_nn_2012_cluster1.ens", h=F)
df.gilman.list[["cluster1a"]] <- read.csv("../gene_lists/gilman_nn_2012_cluster1a.ens", h=F)
df.gilman.list[["cluster1b"]] <- read.csv("../gene_lists/gilman_nn_2012_cluster1b.ens", h=F)
df.gilman.list[["cluster2"]] <- read.csv("../gene_lists/gilman_nn_2012_cluster2.ens", h=F)



### Checking mapping
lapply(df.gilman.list, function(x) sum(!x[,1] %in% df.expression_matrix.clean$ensembl_gene_id)) # unmapped genes
## check that all Gilman clusters are disjoint!
s <- unlist(lapply(df.gilman.list, function(df) as.character(df[,1])))
print(paste("length full list:", length(s))); print(paste("length of unique list:",length(unique(s))))


#### Adding GILMAN to plot #######
df.gilman.melt <- df.expression_matrix.clean.melt
df.gilman.melt$cluster <- NA
### ** OBS - potentially dangerous to use "%in% df.expression_matrix.clean.melt$ensembl_gene_id" - make sure that the df's are on sync (freshly copied) ***
df.gilman.melt[(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.gilman.list[["cluster1a"]][,1]), "cluster"] <- "cluster1a"
df.gilman.melt[(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.gilman.list[["cluster1b"]][,1]), "cluster"] <- "cluster1b"
df.gilman.melt[(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.gilman.list[["cluster2"]][,1]), "cluster"] <- "cluster2"
df.gilman.melt$cluster <- as.factor(df.gilman.melt$cluster)
### Removing NA rows
df.gilman.melt <- df.gilman.melt[complete.cases(df.gilman.melt),]
table(df.gilman.melt$cluster)

str(df.gilman.melt)

########## *** CONSIDER DOING this ddply in the individual scripts
df.gilman.summary.mean <- ddply(df.gilman.melt, c("stage", "cluster"), summarise,
                                mean = mean(value, na.rm=TRUE),
                                sd   = sd(value, na.rm=TRUE))


