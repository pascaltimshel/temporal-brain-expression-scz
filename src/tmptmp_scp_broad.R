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


############################ GGPLOT #########################
p <- ggplot()
#p <- p + 
### plot all lines
p <- p + geom_line(data=df.fit.stage, aes(x=stage, y=t_statistic, group=perm, color=perm))
p
### plot loess smooting
#p + stat_smooth(data=df.fit.stage, aes(x=stage, y=t_statistic, group=1), method="loess",se=T, level=0.95)
### plot sd as ribbon
df.fit.stage.summary <- ddply(df.fit.stage, .(stage), summarise, mean=mean(t_statistic), sd=sd(t_statistic))
p <- p + geom_line(data=df.fit.stage.summary, aes(x=stage, y=mean, group=1), size=1)
p + geom_ribbon(data=df.fit.stage.summary, aes(x=stage, group = 1, ymin=mean-sd, ymax=mean+sd), alpha=0.5, fill="black", color="gray")

### plot summary_stat - DID NOT WORK
p <- ggplot()
p + stat_summary(data=df.fit.stage.summary, aes(x=stage, y=mean, group = 1), fun.data ="mean_cl_normal", geom = "smooth")
p + stat_summary(data=df.fit.stage.summary, aes(x=stage, y=mean, group = 1), fun.data ="mean_sdl", geom = "smooth")

mean_cl_normal(df.fit.stage.summary)

########################### STAT_SUMMARY #######################
#http://stackoverflow.com/questions/13341056/plotting-average-of-multiple-variables-in-time-series-using-ggplot
stat_summary(fun.data = "mean_cl_boot", geom = "smooth")
stat_summary(fun.data ="mean_sdl", mult=1, geom = "smooth") # OBS: mult

#########
#http://docs.ggplot2.org/current/stat_summary.html
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

d + stat_sum_df("mean_cl_boot")
d + stat_sum_df("mean_sdl", mult=1) # OBS: parsing mult using "..."



# Alternatively, you can supply a function that operates on a data.frame.
# A set of useful summary functions is provided from the Hmisc package:
d <- qplot(cyl, mpg, data=mtcars)
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

d + stat_sum_df("mean_cl_normal", geom = "smooth")

