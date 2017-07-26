#!/bin/csh

../Analyze_Metadata.r \
	-f mousePH.metadata.SB.tsv \
	--exclude exclude.list

../Analyze_Metadata.r \
	-f mousePH.metadata.LB.tsv \
	--exclude exclude.list
