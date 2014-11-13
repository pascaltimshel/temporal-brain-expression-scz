
stages = list()
# 1 4-7 pcw Embryonic
# none

# 2A 8-9 pcw Early prenatal
stages[["s2a"]] = c("8 pcw","9 pcw")

# 2B 10-12 pcw Early prenatal
stages[["s2b"]] = c("12 pcw")

# 3A 13-15 pcw Early mid-prenatal
stages[["s3a"]] = c("13 pcw")

# 3B 16-18 pcw Early mid-prenatal
stages[["s3b"]] = c("16 pcw","17 pcw")

# 4 19-24 pcw Late mid-prenatal
stages[["s4"]] = c("19 pcw","21 pcw","24 pcw") 

# 5 25-38 pcw Late prenatal
stages[["s5"]] = c("25 pcw","26 pcw","35 pcw","37 pcw")

# 6 Birth-5 months Early infancy
stages[["s6"]] = c("4 mos")

# 7 6-18 months Late infancy
stages[["s7"]] = c("10 mos")

# 8 19 months-5 yrs Early childhood
stages[["s8"]] = c("1 yrs","2 yrs","3 yrs","4 yrs")

# 9 6-11 yrs Late childhood
stages[["s9"]] = c("8 yrs","11 yrs")

# 10 12-19 yrs Adolescence
stages[["s10"]] = c("13 yrs","15 yrs","18 yrs","19 yrs")

# 11 20-60+ yrs Adulthood
stages[["s11"]] = c("21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")

# RNA-seq expression data
#flag = "microarray"
flag = "rnaseq"
path = paste("/home/data/expression/brainspan/data/141031/",flag,"/",sep="")
exp = read.csv(paste(path,"expression_matrix.csv",sep=""),row.names=1,h=F)
columns = read.csv(paste(path,"columns_metadata.csv",sep=""),h=T,row.names=1)
colnames(exp) = paste(gsub(" ","-",columns$age),columns$structure_acronym,sep="_")
rownames.tmp = unlist(read.csv(paste(path,"rows_metadata.csv",sep=""),h=T,row.names=1)[2])
if (sum(duplicated(rownames.tmp))) { # Remove duplicates genes and genes with no Ensembl IDs
	exp = exp[-which(rownames.tmp==""),]
	rownames.tmp = rownames.tmp[-which(rownames.tmp=="")]
	exp = exp[-which(duplicated(rownames.tmp)),]
	rownames(exp) = rownames.tmp[-which(duplicated(rownames.tmp))]
} else {
	rownames(exp) = rownames.tmp
}

# Function to normalized rows and columns
norm <- function(a) {
	b = a[rowMeans(a)!=0,]	
	return(scale(b-rowMeans(b)/apply(b,1,sd)))
}
exp.norm = norm(exp)

# Construct collapsed expression matrix
count = 0
exp.merged = matrix(ncol=1,nrow=nrow(exp.norm))
header = c()
for(stage in names(stages)) {


	# Find samples for the given stage
	indices = c()
	for (age in stages[[stage]]) {
		indices = c(indices,grep(paste("^",age,"$",sep=""),columns$age,perl=T))
	}
	# .. and structures for these samples
	structures = unique(columns$structure_acronym[indices])
	for (structure in structures) {

		# Merge by structure within stage
		merge_indices = indices[grep(paste("^",structure,"$",sep=""),columns$structure_acronym[indices])]
		print(paste(colnames(exp.norm)[merge_indices]))
		count = count + length(merge_indices)

		# Average expression data
		header = c(header,paste(stage,structure,sep="_"))
		if (length(merge_indices) == 1 ) {
			merged.data = exp.norm[,merge_indices]
		} else {
			merged.data = rowMeans(exp.norm[,merge_indices])
		}
		exp.merged = data.frame(exp.merged,merged.data)
	}
}

# Alle tissues captured?
ncol(exp.norm) == count

exp.merged = exp.merged[,-1]
colnames(exp.merged) = header

print("REMEMBER to add '-' in fron of header")
write.table(exp.merged,paste(path,"brainspan_expression.tab",sep=""),row.names=T,quote=F,sep="\t")



