#!/usr/bin/env python

import glob
import sys
import os

import pdb

################## EXAMPLE file ##################
####### NOTE: THE PART OF THE FILE WE ARE INTERESTED IN IS SORTED!
# CRAP IN BEGINNING OF FILE
# .....
# Locus	LocusDescription	NrGenes	ChrPos	TopGene	TopGeneZScore	TopGenePValue	TopGenePValueCorrectedForLocusSize	TopGeneClosestToSNP	GeneBioType	TopCisEQTLSNP	FDR	NrTPs
# 3	rs8030757	3	chr15-38746328-38857776	ENSG00000259326	2.0192541676372104	0.021730404451074496	0.021730404451074496	false	antisense	-	19.0	-18.0
# 67	rs6459434	1	chr6-9596343-10211841	ENSG00000181355	1.9374177215891561	0.0263471510426946	0.0263471510426946	true	processed_transcript+protein_coding+nonsense_mediated_decay	rs2064082	10.0	-18.0
# 63	rs10168607;rs11688737;rs11892771;rs12712596;rs12712602;rs13019288;rs1368112;rs6749177;rs906806	1	chr2-38658436-38660715	ENSG00000235009	1.7943532000360878	0.036378405861619244	0.036378405861619244	true	lincRNA	-	7.0	-18.0
# 47	rs10756828;rs12683636;rs13289064;rs1339549;rs1416745;rs4436220;rs7868212	1	chr9-16416404-16870841	ENSG00000173068	1.6991712587349996	0.04464346003169906	0.04464346003169906	true	protein_coding+nonsense_mediated_decay+processed_transcript	rs10962474	5.25	-17.0
# ....
# MORE CRAP
###

# Read genes from associated loci
def get_genes(loci_file):
	with open(loci_file, 'r') as infile:
		flag_loci = 0
		locigenes = {}
		for line in infile.readlines():
			words = line.strip("\n").split("\t") 
			#NOTE: if line is empty ('\n'): words will be a list of length 1 (['']). Splitting an empty string with a specified separator returns [''].
	
			if words == [''] and flag_loci: # end of data of interest!
				flag_loci = 0
	
			# Read all loci and their genes
			if len(words) > 1 and words[0] == "Locus" and words[1] == "LocusDescription":
				flag_loci = 1
				continue

			if flag_loci:
				# Gem gener og deres z-score (eller p value)
				TopGene = words[4]
				TopGenePValue = words[6]
				locigenes[TopGene] = TopGenePValue
	return locigenes



n_top_genes = 63
path_in = "/Users/pascaltimshel/p_scz/brainspan/data/schizophrenia_expression_geneprio"
#file_out = "schizophrenia_expression325permutations0to999.genes.prioritized.combined.csv"
file_out = "schizophrenia_expression325permutations0to999.genes.prioritized.top_%s.combined.csv" % n_top_genes

glob_pattern = path_in+"/*.out"
print (glob_pattern)
permfiles = glob.glob(glob_pattern)
permfiles.sort(key=lambda x: int(x.split("perm")[-1].split(".")[0]))
print (len(permfiles))
### EXAMPLE file name
#/Users/pascaltimshel/p_scz/brainspan/schizophrenia_expression_geneprio/schizophrenia_expression_loci_permutation13_325_GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z_affymetrix_geneprioritization.txt

permutation_check_order = []
with open(file_out, 'w') as fout:
	fout.write("permutation,ensembl_gene_id,pval_nominal\n")
	for i, permfile in enumerate(permfiles, start=1):

		## Get permutation number
		#perm_no = permfile.split("permutation")[-1].split("_")[0]
		perm_no = permfile.split("perm")[-1].split(".")[0]
		print "perm_no=%s; permfile=%s" % (perm_no, permfile)
	 	permutation_check_order.append(int(perm_no)) # OBS: CONVERTING TO int 

	 	locigenes = get_genes(permfile)
	 	if len(locigenes) < n_top_genes:
	 		#raise Exception("len(locigenes) < n_top_genes")
	 		 print "ERROR: len(locigenes) < n_top_genes"
	 		 print "len(locigenes)=%s" % len(locigenes)
	 		 raw_input("Press Enter to continue...")

	 	for gene in sorted(locigenes, key=locigenes.get)[:n_top_genes]: #or __getitem__ ---> SORTED BY P-value
	 	#for gene in sorted(locigenes)[:n_top_genes]: #or __getitem__
	 		fout.write("%s,%s,%s\n" % (perm_no, gene, locigenes[gene]))
	 	print "len(locigenes)=%s" % len(locigenes)

if not permutation_check_order == range(0,1000): # gives 0...999
	setdiff1 = set(permutation_check_order).difference(set(range(0,1000)))
	setdiff2 = set(range(0,1000)).difference(set(permutation_check_order))
	setdiff_check = set(range(0,100)).difference(set(range(0,10)))
	print "setdiff1", setdiff1
	print "setdiff2", setdiff2
	print "setdiff_check", setdiff_check
	
	print permutation_check_order

	#print "setdiff1", "|".join(list(setdiff1))
	#print "setdiff2", "|".join(list(setdiff2))
	#print "setdiff_check", "|".join(list(setdiff_check))
	raise Exception("permutation is not in expected order! Fail")
else:
	print "permutations are of correct size and order: 0, 1, ..., 998, 999"

