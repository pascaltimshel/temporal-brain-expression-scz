
####### COOL commands #######
#ggplot_build(<plot>)
#

## expand_limits(x = 0) --> include x=0 in the view of the plot
## scale_y_continuous(expand = c(0,1)) # multiplicative and additive constant used to expand the range of the scales 
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
