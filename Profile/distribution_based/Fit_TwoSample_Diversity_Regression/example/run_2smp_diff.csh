#!/bin/csh

../Fit_PairedDiff_Diversity_Regression.r \
	-s 0108_Stapleton.taxa.genus.no_mcl.summary_table.tsv \
	-p np_sinus.crs_obst.first_visit.NasalCavity.map.tsv \
	-f Metadata.joined.tsv \
	-F SampleID \
	-M covariates.sample_type_analysis \
	-B NP \
	-A SINUS \
	-o npsinus.crs_obst.firstvisit


