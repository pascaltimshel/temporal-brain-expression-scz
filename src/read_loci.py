#!/usr/bin/env python

import glob
import sys
import os


################## Plot histogram of number of loci distributed over permutations (n=1000) ##################
# pascaltimshel@mbp:~/p_scz/brainspan/src$ python read_loci.py | grep 'n_loci=' | perl -lne '/=(\d+)/; print $1' > tmp.tmp_n_loci
### R
# x<-read.table("tmp.tmp_n_loci")
# str(x)
# ggplot(x, aes(x=V1)) + geom_histogram(binwidth=1)
# mean(x$V1)
# range(x$V1)
# > mean(x$V1)
# [1] 316.7
# > range(x$V1)
# [1] 302 325


# Read genes from associated loci
def get_genes(loci_file):
	n_loci = 0
	locigenes = {}
	with open(loci_file,'r') as infile:
		lines = infile.readlines()
		for line in lines:
			words = line.strip().split()
			if words[0] != "#" and words[0] != "snps":
				n_loci += 1
				for g in words[4].split(';'):
					locigenes[g] = g
				if len(words) > 5:
					for g in words[5:len(words)]:
						locigenes[g] = 1
	print "n_loci=%s" % n_loci
	return locigenes # dict
	
path_in = "/Users/pascaltimshel/p_scz/brainspan/schizophrenia_expression325permutations0to999"
file_out = "schizophrenia_expression325permutations0to999.genes.combined.csv"

permfiles = glob.glob(path_in+"/*.txt")
#/Users/pascaltimshel/p_scz/brainspan/schizophrenia_expression325permutations0to999/schizophrenia_expression_loci_permutation0_325.txt
#permutation0_325.txt
#print "\n".join(permfiles)

with open(file_out, 'w') as fout:
	fout.write("permutation,ensembl_gene_id\n")
	for i, permfile in enumerate(permfiles, start=1):
		perm_no = permfile.split("permutation")[-1].split("_")[0]
		print perm_no, permfile
	 	locigenes = get_genes(permfile)
	 	for gene in locigenes:
	 		fout.write("%s,%s\n" % (perm_no, gene))
	 	print len(locigenes)

