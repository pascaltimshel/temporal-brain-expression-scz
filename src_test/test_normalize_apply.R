#df <- data.frame(a=rnorm(50)+5, b=rnorm(50))
#df <- data.frame(a=sample(50), b=sample(50)) #--> standadization of 2 numbers will always sqrt(2)/2
df <- data.frame(a=sample(500), b=sample(500), c=sample(500))
df
apply(df,1, sum)
rowMeans(df)
apply(df,1,sd)
df - rowMeans(df)
df.norm <- (df - rowMeans(df))/apply(df,1,sd)

### check
rowMeans(df.norm)
apply(df.norm,1,sd)

x/c(2,5,30)

apply(df,1,sd)

(410-158.0)/223.445743