xx <- data.frame(group = rep(1:4, 100), a = rnorm(400) , b = rnorm(400) )
head(xx)

func <- function(xx) {
  return(data.frame(COR = cor(xx$a, xx$b)))
}

ddply(xx, .(group), func)
ddply(xx, c("group"), func)
ddply(xx, .(group), summarize,
      corr=cor(a,b)) # USE THIS VERSION


va = rnorm(5)
vb = va*5+rnorm(5)
xx <- data.frame(a=va,b=vb)
xx <- data.frame(a=va,b=vb, c=c(vb[-c(1,2)], 5, NA))
xx
cor(xx, use="complete.obs")
cor(xx, use="pairwise.complete.obs")
cor(xx$a, xx$b)
#"complete.obs" and "pairwise.complete.obs" will give the same


