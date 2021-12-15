#!/bin/csh

../Double_Paired_Difference.r \
	-s ST.summary_table.tsv \
	-m ST.prepost_cat.map.tsv \
	-a pre \
	-b post \
	-n Stool \
	-f ST.paired_diff.questionnaire.metadata.tsv \
\
	-S SAL.summary_table.tsv \
	-M SAL.prepost_cat.map.tsv \
	-A pre \
	-B post \
	-N Saliva \
	-F SAL.paired_diff.questionnaire.metadata.tsv \
\
	-c target_variables \
	-l collection_date_int \
	-w AX \
\
	-o results \
	-d euc
