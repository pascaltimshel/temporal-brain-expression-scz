a <- gl(2, 4, 8)
b <- gl(2, 2, 8, labels = c("ctrl", "treat"))
s <- gl(2, 1, 8, labels = c("M", "F"))
interaction(a, b)
interaction(a, b, s, sep = ":")

?gl

#http://stackoverflow.com/questions/10357768/plotting-lines-and-the-group-aesthetic-in-ggplot2
#http://kohske.wordpress.com/2010/12/27/faq-geom_line-doesnt-draw-lines/

#####################################

#http://docs.ggplot2.org/current/aes_group_order.html

# There are three common cases where the default is not enough, and we
# will consider each one below. In the following examples, we will use a simple
# longitudinal dataset, Oxboys, from the nlme package. It records the heights
# (height) and centered ages (age) of 26 boys (Subject), measured on nine
# occasions (Occasion).

# Multiple groups with one aesthetic
library(nlme)
h <- ggplot(Oxboys, aes(age, height))
# A single line tries to connect all the observations
h + geom_line()

# The group aesthetic maps a different line for each subject
h + geom_line(aes(group = Subject))

# Different groups on different layers
h <- h + geom_line(aes(group = Subject))
# Using the group aesthetic with both geom_line() and geom_smooth()
# groups the data the same way for both layers
h + geom_smooth(aes(group = Subject), method = "lm", se = FALSE)

# Changing the group aesthetic for the smoother layer
# fits a single line of best fit across all boys
h + geom_smooth(aes(group = 1), size = 2, method = "lm", se = FALSE)
h + geom_smooth(size = 2, method = "lm", se = FALSE) # PASCAL TEST: this gives the same as with "group = 1" 
