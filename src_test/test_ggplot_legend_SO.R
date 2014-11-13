#http://stackoverflow.com/questions/16320148/ggplot2-draw-dashed-lines-of-same-colour-as-solid-lines-belonging-to-different-g

x <- c(10, 20, 50, 10, 20, 50)
mean = c(52.4, 98.2, 97.9, 74.1, 98.1, 97.6)
group = c(1, 1, 1, 2,2,2) 
upper = c(13.64, 89, 86.4, 13.64, 89, 86.4)
lower = c(95.4, 99.8, 99.7, 95.4, 99.8, 99.7)
data <- data.frame(x=x,y=mean, group, upper, lower)

ggplot(data, aes(x = x, y= mean, group = as.factor(data$group), colour=as.factor(data$group))) + geom_line() + geom_point() + geom_line(data=data,aes(x=x, y=lower, group = as.factor(data$group), colour=as.factor(data$group), linetype="dotted")) + geom_line(data=data,aes(x=x, y=upper, group = as.factor(data$group), colour=as.factor(data$group), linetype="dotted")) + scale_color_manual(values=c("red", "blue")) +  scale_colour_discrete(name="Groups")

ggplot(data, aes(x = x, y= mean, group = group)) + geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper)) +
  geom_line(aes(y = mean), colour = "Mean")) +
  scale_colour_manual(name = "", values = c("Group1", "Group2"))


p <- ggplot(data, aes(x = x, y= mean, group = as.factor(data$group), colour=as.factor(data$group))) + 
  geom_line() + geom_point() + 
  geom_line(aes(y=lower),linetype="dotted") + 
  geom_line(aes(y=upper),linetype="dotted")+
  scale_color_manual(name="Groups",values=c("red", "blue"))
p
p + guides(colour = guide_legend(override.aes = list(linetype = c(4,2))))



################################
#http://stackoverflow.com/questions/5435883/how-to-merge-colour-and-shape
dataf <- structure(list(Municipality = structure(c(2L, 4L, 10L, 11L, 6L, 8L, 3L, 1L, 5L, 9L, 7L), .Label = c("Boyuibe", "Cabezas", "Camiri", "Charagua", "Cuevo", "Gutierrez", "Huacaya", "Lagunillas", "Machareti", "Vallegrande", "Villa Vaca Guzman"), class = "factor"), Growth = c(3.05, 2.85, 0.14, 1.21, 1.59, 2.35, -0.41, 0.81, 0.9, 2.89, 1.8), Density = c(3.0390920594, 0.260984024187, 5.20069847261, 2.50828556783, 3.43964629267, 3.69768961375, 32.4496626479, 2.06145019368, 4.2139578988, 0.740736713557, 1.67034079825)), .Names = c("Municipality", "Growth", "Density"), class = "data.frame", row.names = c(NA, -11L))

dataf <- dataf[with(dataf, order(Municipality)), ]
# create a new column with values 1 to 6 and same length as Municipality
modulus <- function(x) (x - 1) %% 6 + 1
indeces <- 1:length(dataf$Municipality)
dim(indeces) <- length(dataf$Municipality)
dataf$Shape <- apply(indeces, 1, modulus)
dataf$Shape <- factor(dataf$Shape, levels=unique(dataf$Shape))

plot1 <- ggplot(dataf, aes(x=Density, y=Growth, colour=Municipality,
                           shape=Municipality))
plot1 <- plot1 + geom_point(size=3)
plot1
#plot1 <- plot1 + scale_colour_discrete() 
plot1
plot1 <- plot1 + scale_shape_manual(values=as.numeric(dataf$Shape))
plot1


############################## PLAYING WITH GGPLOT LEGENDS #######################
g <- guide_legend("Legend")
p + guides(colour=g, linetype=g)

p + guides(color=guide_legend(override.aes=list(color=c("red","blue", "black"),linetype=c(1,0,1))))

p + labs(colour = "Cylinders") # The labs function also modifies legend labels


guides(fill=guide_legend(title=NULL)) # No legend

guides(colour = guide_legend(override.aes = list(size=3)))


#http://stackoverflow.com/questions/18060116/adding-legend-to-ggplot-when-lines-were-added-manually
geom_line(aes(x, y, color="My Line"), data=line.data) +
  scale_color_manual(values=c("setosa"="blue4", "versicolor"="red4",
                              "virginica"="purple4", "My Line"="gray"))

