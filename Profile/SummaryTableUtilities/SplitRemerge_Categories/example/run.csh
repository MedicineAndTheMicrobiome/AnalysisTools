#!/bin/csh

../SplitRemerge_Categories.r \
	#-i example_EC.summary_table.tsv \
	-i example_test.summary_table.tsv \
	-o example_EC.splitmerged.summary_table.tsv
