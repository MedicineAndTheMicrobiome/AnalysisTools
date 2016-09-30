#!/bin/csh

./Plot_Longitudinal_Samples.r \
	-i mcdyer.taxa.genus.summary_table.tsv \
	-t mcdyer.time_info.tsv \
	-o /media/sf_SSDScratch/mcdyer.taxa

./Plot_Longitudinal_Samples.r \
	-i mcdyer.otu.97.summary_table.tsv \
	-t mcdyer.time_info.tsv \
	-o /media/sf_SSDScratch/mcdyer.otu
