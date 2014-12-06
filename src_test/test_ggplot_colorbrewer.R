library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer")
library(RColorBrewer)

#### LIST OF R colors!
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
#

#http://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html#rcolorbrewer
display.brewer.all()
#sequential: great for low-to-high things where one extreme is exciting and the other is boring, like (transformations of) p-values and correlations (caveat: here I'm assuming the only exciting correlations you're likely to see are positive, i.e. near 1)
#qualitative: great for non-ordered categorical things -- such as your typical factor, like country or continent. Note the special case "Paired" palette; example where that's useful: a non-experimental factor (e.g. type of wheat) and a binary experimental factor (e.g. untreated vs. treated).
#diverging: great for things that range from "extreme and negative" to "extreme and positive", going through "non extreme and boring" along the way, such as t statistics and z scores and signed correlations

display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = "Dark2")

rm(list=ls())

###################################
# Two variables
df <- read.table(header=T, text='
 cond yval
    A 2
    B 2.5
    C 1.6
')

df2 <- read.table(header=T, text='
 cond1 cond2 yval
    A      I 2
    A      J 2.5
    A      K 1.6
    B      I 2.2
    B      J 2.4
    B      K 1.2
    C      I 1.7
    C      J 2.3
    C      K 1.9
')


cmap <- c("A"="gray", "B"="#d7191c", "C"="black")
p <- ggplot(df, aes(x=cond, y=yval, fill=cond)) + geom_bar(stat="identity") 
p
p + scale_fill_manual(name="MyTest", values=cmap)
p1 <- p + scale_fill_brewer(palette="Set1")
p1 
p1 + scale_fill_manual(name="MyTest", values=cmap) # ---> gives the same as p.
### Conclusion: this uses the 




#http://sape.inf.usi.ch/quick-reference/ggplot2/colour
#http://zevross.com/blog/2014/08/04/beautiful-plotting-in-r-a-ggplot2-cheatsheet-3/#working-with-colors

scale_color_manual(values=c("dodgerblue4", "darkolivegreen4",
                            "darkorchid3", "goldenrod1"))
#scale_color_brewer(palette="Set1")