#!/bin/csh

../Chemical_Formula_Search.r \
	-i example.input.tsv \
	-c Formula \
	#-d ../species.txt \
	-m ../species.txt.mat \
	#-d ../species.test.txt \
	#-d species.db_subset.txt \
	#-m species.db_subset.txt.mat \
	#-m ../species.test.txt.mat \
	-o results
