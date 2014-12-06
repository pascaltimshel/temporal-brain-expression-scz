#!/usr/bin/env python2.7

import os,pdb

# Function to generate HGNC to ensembl mapping
def get_mapping(flag='hgnc'):
	"""Default is to map ENSEML to HGNC"""
	e2h = {}
	h2e = {}
	mapping_file = open(os.path.join('/Users/pascaltimshel/p_scz/brainspan/data/gene_mapping/mapping_ENSG_HGNC_Entrez_v73.tab') , 'r')
	lines = mapping_file.readlines()[1:]
	for line in lines:
		words = line.strip().split('\t')
		if "ENS" in words[0]:
			if flag == 'hgnc':
				if len(words) > 2 and words[2] != '':
					e2h[words[0]] = words[2]
					h2e[words[2]] = words[0]
			elif flag == 'entrez':
				if len(words) > 1 and words[1] != '':
					e2h[words[0]] = words[1]
					h2e[words[1]] = words[0]
	mapping_file.close()
	return e2h, h2e

