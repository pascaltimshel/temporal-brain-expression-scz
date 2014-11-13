#!/usr/bin/env python

with open("variable.txt") as f:
	lines = f.readlines()

for line in lines:
	line = line.strip("\n")
	fields = line.split("_")
	if not len(fields) == 3: print "error; %s" % line

