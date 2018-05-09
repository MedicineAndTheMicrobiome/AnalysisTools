#!/bin/csh

Fit_Diversity_Regression.r \
	-s sample_data/allSamples.padded.genus.summary_table.xls \
	-f sample_data/factorfile.patient.txt \
	-r sample_data/reference_levels.txt
