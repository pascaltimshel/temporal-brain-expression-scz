library(biomaRt)
x.marts <- listMarts()
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
x.att <- listAttributes(ensembl) # start_position, end_position
x.filt <- listFilters(ensembl) # ensembl_gene_id
filterType("ensembl_gene_id",ensembl)

#test_results <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position"), filters = "ensembl_gene_id", values = as.character(df.gilman.list[["cluster1"]][,1]), mart = ensembl)

map.ensembl <- function(list.name) {
  cat(list.name,"\n")
  gene_list <- list2run[[list.name]][,1]
  df <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position"), filters = "ensembl_gene_id", values = gene_list, mart = ensembl)
  df$gene_length <- df$end_position - df$start_position
  n_unmapped <- setdiff(gene_list, df$ensembl_gene_id) # setdiff(x,y) --> genes in x but NOT in y
  cat("unmapped genes: ", length(n_unmapped), "\n", sep="")
  return(df)
}

### Run function
#list2run <- df.gilman.list
list2run <- list(BrainSpan=as.data.frame(df.expression_matrix.clean$ensembl_gene_id))
results <- llply(names(list2run), map.ensembl)
names(results) <- names(list2run)
names(results)

################ Prossesing results ####################
df.results <- results[[1]]
### initialyzing data frame with NA
df.clean.gene_length.insync <- data.frame(ensembl_gene_id=df.expression_matrix.clean$ensembl_gene_id, gene_length=NA)
sum(!df.clean.gene_length.insync$ensembl_gene_id %in% df.results$ensembl_gene_id) # unmapped
### ** takes a long time ** - important to loop over df.results$ensembl_gene_id. In this way we know that all genes are present in df.clean.gene_length.insync as well
for (x in df.results$ensembl_gene_id) {df.clean.gene_length.insync[df.clean.gene_length.insync$ensembl_gene_id==x,"gene_length"] <- df.results[df.results$ensembl_gene_id==x, "gene_length"]}
sum(is.na(df.clean.gene_length.insync$gene_length))



### test - run this a couple of times
testint <- sample.int(15000, 1)
testid <- df.clean.gene_length.insync[testint, "ensembl_gene_id"]
df.clean.gene_length.insync[df.clean.gene_length.insync$ensembl_gene_id==testid,]
df.results[df.results$ensembl_gene_id==testid, ]


#### WRITE OUT RESULTs
write.csv(df.clean.gene_length.insync, "df.clean.gene_length.insync.csv", row.names=F)


