#!/bin/csh

../Export_ALR_Values.r \
	-s mousePH.taxa.genus.no_mc.rep_mrg.gt750.exp_only.summary_table.tsv \
	-p 35 \
	-o mousePH.taxa.genus.no_mc.rep_mrg.gt750.exp_only \
	-x ";"
