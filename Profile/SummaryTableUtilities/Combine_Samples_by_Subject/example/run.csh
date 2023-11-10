#!/bin/csh

../Combine_Samples_by_Subject.r \
	-i 0180.SEEDS.16S.DepETA.summary_table.tsv \
	-m 0180.SEEDS.16S.DepETA.alr.top15.visit_num.map.tsv \
	-o combined.v1v2 \
	-l v1,v2

../Combine_Samples_by_Subject.r \
	-i 0180.SEEDS.16S.DepETA.summary_table.tsv \
	-m 0180.SEEDS.16S.DepETA.alr.top15.visit_num.map.tsv \
	-o combined.all 

