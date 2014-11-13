Reduce(intersect,  lapply(df.gilman.list, function(df) as.character(df[,1])))

tmp.list = list(data.frame(v1=letters[1:5]), data.frame(v1=letters[6:10]), data.frame(v1=letters[c(2,7)]))
Reduce(intersect, lapply(tmp.list, function(df) as.character(df[,1])))

lapply(tmp.list, function(df) as.character(df[,1]))

if (length(unique(s)) != length(s)) {print("STOP"); print(paste("length full list:", length(s))); print(paste("length of unique list:",length(unique(s))))} 