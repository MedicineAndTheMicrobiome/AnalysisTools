#!/bin/csh

../calculate_new_columns.r \
	-i grads.subset.tsv \
	-f formulas \
	-o grads.subset.extra.tsv
	

