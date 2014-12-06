
####### COOL commands #######
#ggplot_build(<plot>)
#
attributes(p)
p$coordinates #p$coordinates$limits
attributes(p$scales)

##################### BUILD ##############
ggplot_build(p)
str(ggplot_build(p))
names(ggplot_build(p)) 
#1) "data"--> a list of data frames (one for each layer)
#2) "panel"-->panel object, which contain all information about axis limits, breaks etc.
#3) "plot"--> builds plot
?ggplot_build
########### Useful for getting axis limits. Use in combination with coord_cartesian(ylim=c(0, 150))
ggplot_build(p)$panel$ranges[[1]]$x.range
ggplot_build(p)$panel$ranges[[1]]$y.range
##########

## expand_limits(x = 0) --> include x=0 in the view of the plot
## scale_y_continuous(expand = c(0,1)) # multiplicative and additive constant used to expand the range of the scales 
# ^^^ ----> DEFAULT: scale_y_continuous(expand = c(0.05, 0)) #c(0.05, 0)
##^^^ ----> consider using 
p.x.range <- ggplot_build(p)$panel$ranges[[1]]$x.range
p.y.range <- ggplot_build(p)$panel$ranges[[1]]$y.range
p <- p + coord_cartesian(xlim=p.x.range*1.05, ylim=p.y.range*1.05)


## scale_y_continuous(limits = c(0, NA)) # same as adding "+ ylim(c(0, NA))"
  # ---> this will EXCLUDE observations from the data set
## coord_cartesian(ylim=c(0, 150)) # HARD cut-off. 

### INFO:
# http://docs.ggplot2.org/current/continuous_scale.html

## segment info
# http://stackoverflow.com/questions/9085104/is-there-a-way-to-limit-vline-lengths-in-ggplot2


########### EXAMPLE data ##########

tmp.null_distribution.t_statistic <- rnorm(1000, 3, 1)

########### EXAMPLE COMMANDS #########
q <- qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2)
q
attributes(q)
q$coordinates$limits
q$coordinates

attributes(q$scales)
ggplot_build(q)

qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + scale_y_continuous(expand = c(0,1)) # multiplicative and additive constant used to expand the range of the scales 
qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + coord_cartesian(ylim=c(0, 150)) # HARD cut-off. 

qplot(tmp.null_distribution.t_statistic) + geom_vline(xintercept=c(tmp.obs.t_statistic), linetype="dotted", size=2) + scale_y_continuous(limits = c(0, NA)) # same as adding "+ ylim(c(0, NA))"

##### using geom_segment() - DOES NOT solve issue
qplot(tmp.null_distribution.t_statistic) + geom_segment(aes(x = tmp.obs.t_statistic, xend=tmp.obs.t_statistic, y = 0, yend = Inf), linetype="dotted", size=2) 

#############################

pshare <- data.frame()

for (i in 1:365) {
  pshare <- rbind(pshare,c(i, pbirthday(i,365,coincident=3)))
}

names(pshare) <- c("number","probability")

probs <- c(0.25, 0.50, 0.75)
marks <- data.frame(probability = probs,
                    number = sapply(probs, qbirthday, classes=365, coincident=3))

qplot(number,probability,data=subset(pshare,probability<0.99)) +
  geom_segment(data=marks, aes(xend=-Inf, yend=probability)) +
  geom_segment(data=marks, aes(xend=number, yend=-Inf))
