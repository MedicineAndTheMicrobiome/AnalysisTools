#!/bin/csh

~/git/AnalysisTools/Profile/distance_based/Distance_PML/Distance_PML.r \
	-d mousePH.taxa.genus.sb.man.distmat \
	-f mousePH.metadata.SB.tsv \
	-o mousePH.taxa.genus.SB.man \
	--include include_short.txt  | tee stderr.man.sb 

