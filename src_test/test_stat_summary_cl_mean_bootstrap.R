
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

