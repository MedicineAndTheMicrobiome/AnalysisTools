#!/bin/csh

../Cluster_Transistion_TwoProfiles.r \
	-s SAL_PrePost_COVID.summary_table.tsv \
	-m SAL.PrePost_cat.map.tsv \
	-A pre \
	-B post \
	-o SAL \
	-f SAL.metadata.tsv
