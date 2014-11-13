library(ggplot2)
library(plyr) # must be explicitly loaded

DF <- data.frame(a = letters[1:10], b = 1:10)
ggplot(DF, aes(x = b, y = b)) + 
  geom_point(subset = .(a == 'a'), colour = 'blue') + 
  geom_point(subset = .(a == 'c'), colour = 'green') +
  geom_line()

.(a == 'c')
