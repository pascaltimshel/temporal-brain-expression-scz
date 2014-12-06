load("null_RData_broad_marray_priority.RData")
x<-list.par_analysis[[1]]
x<-list.par_analysis[[1000]]
x$df.null.mapping


x$list.null.fits$list.fit.stage[1]


############################################################

############################# EXTRACT BROAD DATA #################################
############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal <- lapply(list.null.fits, "[[", "fit.natal")
list.null.fit.prioritized.higher <- lapply(list.null.fits, "[[", "fit.prioritized.higher")
list.null.list.fit.stage <- lapply(list.null.fits, "[[", "list.fit.stage")
############## Extracting from list.null.list.fit.stage
df.fit.stage<-ldply(list.null.list.fit.stage, function(y) {ldply(y, function(x) {x$statistic})})
list.null.list.fit.stage[[1]]
names(list.null.list.fit.stage)

#### Investigation of how ldply "gets it right"
class(list.null.list.fit.stage[[1]][["s2a"]]$statistic)
attributes(list.null.list.fit.stage[[1]][["s2a"]]$statistic)
as.data.frame(list.null.list.fit.stage[[1]][["s2a"]]$statistic)
names(list.null.list.fit.stage[[1]])
attributes(list.null.list.fit.stage[[1]]) # it contains the "$split_labels" attribute. This must be how ldply can re-combine the proper names automatically (stage ID)
ldply(list.null.list.fit.stage[[1]], function(fit.stage) {fit.stage$statistic}, .id="bla")
# the BELOW WORKS!
df.fit.stage<-ldply(seq_along(list.null.list.fit.stage), function(i) {ldply(list.null.list.fit.stage[[i]], function(fit.stage) {data.frame(t_statistic=fit.stage$statistic, perm=i)})})

############## Combining summary data frames
df.null.mean.summary <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary <- ldply(list.null.median.summary) # COMBINING list of data frames

