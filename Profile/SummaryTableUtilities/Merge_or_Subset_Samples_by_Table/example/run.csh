#!/bin/csh

../Merge_or_Subset_Samples_by_Table.r \
	-i example.summary_table.tsv \
	-m example.collapse.map.tsv \
	-M 3000 \
	-o example
