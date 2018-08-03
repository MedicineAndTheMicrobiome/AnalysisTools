#!/bin/csh

../SubsetSamples_by_FactorInfo.r \
	-f grads.subset.extra.tsv \
	-c 'bpHorL=="High" & demog.gender=="Male"' \
	-o grads.subset.extra.HighBP_Male.tsv
	

