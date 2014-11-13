#!/usr/bin/env python

import glob
import sys
import os


################## EXAMPLE file ##################
# Locus	Nr of genes in locus	Chromosome and position	Ensembl Gene ID	Gene symbol	Nominal P value	Gene closest to lead SNP	Gene bio-type	Top cis eQTL SNP (Westra et al. Nature Genetics 2014)	False discovery rate < 5%
# rs7333291	1	chr-1:54689924-54707001	ENSG00000234787	13	0.01017358	TRUE	LINC00458	-	No
# rs2702835	2	chr-1:6693076-6742798	ENSG00000245857	8	0.013830601	TRUE		-	No
# rs1013040	1	chr1:12995477-13334272	ENSG00000237515	16	0.014198556	TRUE	SHISA9	-	No
# rs1309821	2	chr-1:64961755-65167553	ENSG00000123213	5	0.018335378	TRUE	NLN	rs2161610	No
# rs2652572	1	chr1:17684693-17955966	ENSG00000249937	5	0.018707446	TRUE		-	No
# rs1540560	1	chr-1:132284871-133402414	ENSG00000183715	11	0.020568553	TRUE	OPCML	rs2509236	No


# Read genes from associated loci
def get_genes(loci_file):
	locigenes = {}
	with open(loci_file,'r') as infile:
		lines = infile.readlines()[1:] # SKIPPING HEADER
		for line in lines:
			words = line.strip().split()
			gene = words[3]
			pval_nominal = float(words[5])
			locigenes[gene] = pval_nominal
	return locigenes # dict
	
path_in = "/Users/pascaltimshel/p_scz/brainspan/schizophrenia_expression_geneprio"
#file_out = "schizophrenia_expression325permutations0to999.genes.prioritized.combined.csv"
file_out = "schizophrenia_expression325permutations0to999.genes.prioritized.top54.combined.csv"

permfiles = glob.glob(path_in+"/*.txt")
### EXAMPLE file name
#/Users/pascaltimshel/p_scz/brainspan/schizophrenia_expression_geneprio/schizophrenia_expression_loci_permutation13_325_GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z_affymetrix_geneprioritization.txt

with open(file_out, 'w') as fout:
	fout.write("permutation,ensembl_gene_id,pval_nominal\n")
	for i, permfile in enumerate(permfiles, start=1):
		perm_no = permfile.split("permutation")[-1].split("_")[0]
		print perm_no, permfile
	 	locigenes = get_genes(permfile)
	 	#for gene in sorted(locigenes, key=locigenes.get): #or __getitem__
	 	for gene in sorted(locigenes, key=locigenes.get)[:54]: #or __getitem__
	 		fout.write("%s,%s,%s\n" % (perm_no, gene, locigenes[gene]))
	 	print len(locigenes)

