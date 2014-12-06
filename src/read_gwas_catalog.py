#!/usr/bin/env python

import glob
import sys
import os

# custom modules
#import mapping_ens_hgnc

import pdb

################## EXAMPLE file ##################
# Locus	Nr of genes in locus	Chromosome and position	Ensembl Gene ID	Gene symbol	Nominal P value	Gene closest to lead SNP	Gene bio-type	Top cis eQTL SNP (Westra et al. Nature Genetics 2014)	False discovery rate < 5%
# rs1108842;rs2590838	11	chr3:52529354-52865495	ENSG00000010327	STAB1	0.005508607708539853	false	protein_coding+retained_intron	rs4687545	No
# rs1108842;rs2590838	11	chr3:52529354-52865495	ENSG00000055955	ITIH4	0.0068625116561194435	false	protein_coding+retained_intron+processed_transcript	rs17331178	No
###

# Column_no
# 0: Locus
# 1: Nr of genes in locus
# 2: Chromosome and position
# 3: Ensembl Gene ID
# 4: Gene symbol
# 5: Nominal P value
# 6: Gene closest to lead SNP
# 7: Gene bio-type
# 8: Top cis eQTL SNP (Westra et al. Nature Genetics 2014)
# 9: False discovery rate < 5%

# Read genes from associated loci
def get_genes(loci_file):
	with open(loci_file, 'r') as infile:
		flag_loci = 0
		locigenes = {}
		for line in infile.readlines()[1:]: # SKIPPING HEADER
			words = line.strip("\n").split("\t") 
			#NOTE: if line is empty ('\n'): words will be a list of length 1 (['']). Splitting an empty string with a specified separator returns [''].
			gene_id = words[3]
			FDR_pass = words[9]
			if FDR_pass == "Yes":
				locigenes[gene_id] = "prioritized"
			elif FDR_pass == "No":
				pass
				#locigenes[gene_id] = "associated" # INCLUDE ASSOCIATED GENES!
			else:
				raise Exception("Got unexpected value for FDR_pass. Expect 'Yes' or 'No'.")
	return locigenes





path_in = "/Users/pascaltimshel/p_scz/brainspan/data/Old Cofunc GWAS Catalog"
file_out = "GWAS_catalog_prioritized_genes.csv"
#file_out = "GWAS_catalog_associated_genes.csv"

glob_pattern = path_in+"/*"
print (glob_pattern)
permfiles = glob.glob(glob_pattern)
permfiles.sort()
print "len(permfiles)=%s" % len(permfiles)

### EXAMPLE file name
#/Users/pascaltimshel/p_scz/brainspan/data/Old Cofunc GWAS Catalog/Adiponectin_levels_affymetrix_GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary_geneprioritization.txt

phenotypes_passing_criteria = []
with open(file_out, 'w') as fout:
	fout.write("permutation,ensembl_gene_id,pval_nominal,n_genes\n")
	for i, permfile in enumerate(permfiles, start=1):

		## Get PHENOTYPE NAME
		perm_no = os.path.basename(permfile).split("affymetrix")[0].rstrip("_") # removing trailing "_"
		print "perm_no=%s; permfile=%s" % (perm_no, permfile)

	 	locigenes = get_genes(permfile)
	 	print "len(locigenes)=%s" % len(locigenes)
	 	if len(locigenes) < 2: # you need at least 2 genes to do the t-test: filter out 
	 		print "did not pass filter of at least one (prioritized/associated) gene in loci: skipping..."
	 		continue

	 	phenotypes_passing_criteria.append(perm_no)
	 	for gene in sorted(locigenes):
	 		fout.write("%s,%s,%s,%s\n" % (perm_no, gene, locigenes[gene], len(locigenes)))
	 	
print "phenotypes_passing_criteria=%s" % len(phenotypes_passing_criteria)
print "\n".join(phenotypes_passing_criteria)


