#install.packages("VennDiagram")
library(VennDiagram)

####################################################################
grid.newpage()
venn.plot <- draw.pairwise.venn(100, 70, 30, c("First", "Second"))
grid.draw(venn.plot)


###################################################################
grid.newpage()
venn.plot <- draw.pairwise.venn(area1        = 100,
                                area2        = 70,
                                cross.area   = 68,
                                scaled       = F,
                                category     = c("First", "Second"),
                                fill         = c("blue", "red"),
                                alpha        = 0.3,
                                lty          = "blank",
                                cex          = 2,
                                cat.cex      = 2,
                                cat.pos      = c(285, 105),
                                cat.dist     = 0.09,
                                cat.just     = list(c(-1, -1), c(1, 1)),
                                ext.pos      = 30,
                                ext.dist     = -0.05,
                                ext.length   = 0.85,
                                ext.line.lwd = 2,
                                ext.line.lty = "dashed")
grid.draw(venn.plot)

######################### Three-way Venn #####################################
grid.newpage()
venn.plot <- draw.triple.venn(65, 75, 85, 35, 15, 25, 5, c("First", "Second", "Third"))
grid.draw(venn.plot)
######################### Three-way Venn - pretty #####################################
grid.newpage()
venn.plot <- draw.triple.venn(area1    = 65,
                              area2    = 75,
                              area3    = 85,
                              n12      = 35,
                              n23      = 15,
                              n13      = 25,
                              n123     = 5,
                              category = c("First", "Second", "Third"),
                              cat.pos  = c(0, 40, 250),
                              cat.dist = c(0.05, 0.05, 0.05),
                              fill     = c("blue", "red", "green"),
                              alpha    = 0.3,
                              lty      = "blank",
                              cex      = 2,
                              cat.cex  = 2,
                              cat.col  = c("blue", "red", "green"))
grid.draw(venn.plot)


grid.draw(draw.triple.venn(1,2,3,0,0,0,0))


######################### Venneuler #####################################
#install.packages("venneuler")
library(venneuler)

### ---> NO TEXT ANNOTATION

vd <- venneuler(c(A=0.3, B=0.3, C=1.1, "A&B"=0.1, "A&C"=0.2, "B&C"=0.1 ,"A&B&C"=0.1))
plot(vd)
# same as c(A=1, `A&B&C`=1, C=1)
m <- data.frame(elements=c("1","2","2","2","3"), sets=c("A","A","B","C","C"))
v <- venneuler(m)
plot(v)
m <- as.matrix(data.frame(A=c(1.5, 0.2, 0.4, 0, 0),
                          B=c(0 , 0.2, 0 , 1, 0),
                          C=c(0 , 0 , 0.3, 0, 1)))
# without weights
v <- venneuler(m > 0)
plot(v)
