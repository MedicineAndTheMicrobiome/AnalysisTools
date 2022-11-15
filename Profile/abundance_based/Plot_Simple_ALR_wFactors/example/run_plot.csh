#!/bin/csh

../Plot_Simple_ALR_wFactors.r \
	-i genus.summary_table.tsv \
	-f factors.tsv \
	-M var_subset.lst \
	-x ";"

