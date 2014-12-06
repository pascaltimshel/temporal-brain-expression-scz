#!/usr/bin/env python

import glob
import sys
import os

# custom modules
import mapping_ens_hgnc

import pdb



def run_mapping():
	""" Map ENSEMBL genes to HGNC. 
	The output order will be the same as the input order.
	Genes not found in the mapping dict will get the empty string ("") value.
	NOTE THAT THE HEADER IS BEING SKIPPED!
	"""
	with open(file_in, 'r') as f_in, open(file_out, 'w') as f_out:
		print "NOTE: skipping header of file_in=%s" % file_in
		lines = f_in.readlines()[1:]
		print "Got %s genes to map.." % len(lines)
		for line in lines:
			words = line.strip().split('\t')
			gene_ens = words[0]
			if gene_ens not in e2h:
				print "Could not find gene_ens %s in e2h" % gene_ens
			gene_hgnc = e2h.get(gene_ens, "")
			f_out.write( "%s,%s\n" % (gene_ens, gene_hgnc) )

###################################### SCZ prioritized genes ######################################
file_in = "/Users/pascaltimshel/p_scz/brainspan/gene_lists/gene_prioritization.txt"
file_out = "/Users/pascaltimshel/p_scz/brainspan/gene_lists/gene_prioritization.csv"
(e2h, h2e) = mapping_ens_hgnc.get_mapping(flag='hgnc') # make the mapping dicts in global variable scope
run_mapping()



