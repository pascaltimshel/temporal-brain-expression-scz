#! /usr/bin/python 

# Mapping
meta_data = {}
with open("/home/data/expression/brainspan/data/141031/microarray/columns_metadata.csv",'r') as infile:
	for line in infile.readlines()[1:]:
		words = [x.replace("\"","") for x in line.strip().split(',')]
		meta_data["%s"%(words[6])] = "%s"%(words[7])

with open("/home/projects/depict/results/schizophrenia/brainspan/microarray/daner_PGC_SCZ52_0513a.gz.p4.clump.areator.sorted.1mhc1_1kg_tissues_r2Threshold0.5-0cMPerMb-0KbExtra.tab",'r') as infile:
	for line in infile.readlines()[1:]:
		words = line.strip().split('\t')
		st = words[1].split("_")[1]
		print("%s\t%s"%(line.strip(),meta_data[st]))
