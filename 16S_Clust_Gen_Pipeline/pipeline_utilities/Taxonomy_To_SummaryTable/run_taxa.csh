#!/bin/csh

Convert_MothurTaxonomy_to_SummaryTable.pl \
	-t sample_list.unique.good.filter.unique.precluster.pick.reference.wang.taxonomy \
	-n sample_list.unique.good.filter.unique.precluster.pick.names \
	-g sample_list.good.pick.groups \
	-o sample_list
