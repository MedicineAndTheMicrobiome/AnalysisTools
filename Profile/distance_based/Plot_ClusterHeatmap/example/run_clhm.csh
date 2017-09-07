#!/bin/csh

../Plot_ClusterHeatmap.r \
	-i grads.genus.balkit.summary_table.tsv \
	-f grads.genus.balkit.noNAs.tsv \
	-d man \
	-M clean.variables
