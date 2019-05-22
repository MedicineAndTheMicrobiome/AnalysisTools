#!/bin/csh

../Annotate_OTU_SummaryTable_Genera.pl \
	-m sample_list.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy \
	-s sample_list.otu.97.summary_table.tsv \
	-o sample_list.otu.97.genus.summary_table.tsv
