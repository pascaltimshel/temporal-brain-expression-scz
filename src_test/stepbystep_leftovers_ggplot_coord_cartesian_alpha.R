###################################### STEP 1 - THIS WORK!!! ################################
p <- ggplot()
### Adding mean Prioritized
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1), alpha=0)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),alpha=0)
### Adding mean ALL (df.all.sem)
p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1), alpha=0)
p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),alpha=0)
p
######### Adding x-tickmarks for stage + more
p <- do_stage_converter(p)
p

###################################### COMMANDS ################################
str(ggplot_build(p))
names(ggplot_build(p)) 
#ggplot_build(p)$panel
#ggplot_build(p)$panel$x_scales

### Save scales
p.x.range <- ggplot_build(p)$panel$ranges[[1]]$x.range
p.y.range <- ggplot_build(p)$panel$ranges[[1]]$y.range

###################################### STEP 1 - This gave the wrong limits. The x-axis was also distorted! ################################
p <- ggplot(data=df.all.sem, aes(x=stage, y=mean1, group=1))
p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color='DUMMY'), alpha=0) # THIS ALSO WORKS.
#p <- p + coord_cartesian(xlim=p.x.range*1.05, ylim=p.y.range*1.05)
#p <- p + coord_cartesian(xlim=p.x.range+0.05, ylim=p.y.range+0.05)
p
######### Adding x-tickmarks for stage + more
p <- do_stage_converter(p)
p

