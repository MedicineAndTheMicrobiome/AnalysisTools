#!/bin/csh

../TimeCorrected_PairedDifference.r \
	-s ST.summary_table.tsv \
	-m ST.prepost_cat.map.tsv \
	-a pre \
	-b post \
	-f ST.paired_diff.questionnaire.metadata.tsv \
	-c target_variables \
	-l collection_date_int \
	-o ST \
	-w A \
	-d man

if(1)then

../TimeCorrected_PairedDifference.r \
	-s SAL.summary_table.tsv \
	-m SAL.prepost_cat.map.tsv \
	-a pre \
	-b post \
	-f SAL.paired_diff.questionnaire.metadata.tsv \
	-c target_variables \
	-l collection_date_int \
	-o SAL \
	-w A \
	-d man

endif
