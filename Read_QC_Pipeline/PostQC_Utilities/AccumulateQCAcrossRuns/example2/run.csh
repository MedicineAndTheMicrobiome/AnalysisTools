#!/bin/csh

../AccumulateQCAcrossRuns.r \
	#-i target_list.tsv \
	-i target_list_1c.tsv \
	-o output.matrix
