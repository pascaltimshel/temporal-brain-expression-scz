
######## ddply, t-test
f <- function(d) t.test(value ~ treatment, data = d)$p.val
ddply(DF, .(gender, variable), f)


########################### Random file generation ####################
path.out <- "/Users/pascaltimshel/p_scz/brainspan/null_gene_list"
set.seed(1)
n_perm = 1000
for (i in 1:n_perm) {print(i)}
sample(df.expression_matrix.clean$ensembl_gene_id, 5)


######### BUGS ###########
max(df.expression_matrix.clean.melt$value)
which.max(df.expression_matrix.clean.melt$value)
which(df.expression_matrix.clean.melt > 4.493551e+264, arr.ind = TRUE)
?which
df.expression_matrix.clean.melt[1336279,]