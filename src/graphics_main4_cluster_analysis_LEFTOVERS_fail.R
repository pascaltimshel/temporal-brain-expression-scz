

##################################################################


plot_gene <- function(df) {
  df.sub <- df
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
  #p <- p + facet_grid(cluster ~ facet_label) # rows ~ columns
  return(p)
}

plot_cluster <- function(df, cluster_index) {
  list_of_cluster_plots <- list()
  for (gene in unique(df[df$cluster==cluster_index, "ensembl_gene_id"])) {
    df.tmp.sub <- df[as.character(df$ensembl_gene_id)==gene, ]
    p <- plot_gene(df.tmp.sub)
    print(class(p))
    #list_of_cluster_plots <- c(list_of_cluster_plots, gene=p)
    list_of_cluster_plots[[gene]] <- p
    #print(class(list_of_cluster_plots))
  }
  #return(p)
  return(list_of_cluster_plots)
}

unique(df.sub$cluster)
list_of_cluster_plots <- plot_cluster(df.sub, 2)
class(list_of_cluster_plots)
length(list_of_cluster_plots)
attributes(list_of_cluster_plots)


list_of_cluster_plots[1:2]

multiplot(plotlist=list_of_cluster_plots)
multiplot(plotlist=list_of_cluster_plots[1:2])


