#!/bin/csh

../Rename_SummaryTable_Categories_wCollapse.r \
	-i example.summary_table.tsv \
	-m model_cat.map \
	-o example.renamed.summary_table.tsv
