library(reshape2)
tips

sp <- ggplot(tips, aes(x=total_bill, y=tip/total_bill)) + geom_point(shape=1)
sp

######### Make label function
mf_labeller <- function(var, value){
  cat(class(value), "\n") #---> factor
  cat(class(var), "\n") # ---> character
  value <- as.character(value)
  cat(value, "\n")
  cat(var, "\n")
  
  if (var=="sex") { 
    value[value=="Female"] <- "Woman"
    value[value=="Male"]   <- "Man"
  }
  return(value)
}

sp + facet_grid(. ~ sex, labeller=mf_labeller)



#######################################

######### Make label function
mf_labeller <- function(var, value){
  value <- as.character(value) # transforming factor into character
  if (var=="ensembl_gene_id") {
    for (value_i in value) {
      row_idx <- match(value_i, df.map.prioritized[,1])
      hgnc_symbol <- as.character(df.map.prioritized[row_idx,2])
      replacement_string <- paste(value_i, "=", hgnc_symbol, sep="")
      value[value==value_i] <- replacement_string
    }
  }
  return(value)
}

#sp + facet_grid(. ~ sex, labeller=mf_labeller)

