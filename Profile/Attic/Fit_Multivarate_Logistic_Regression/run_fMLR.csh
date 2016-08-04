#!/bin/csh

Fit_Multivarate_Logistic_Regression.r \
	-s sample_data/allSamples.padded.genus.summary_table.xls \
	-f sample_data/factorfile.patient.txt 
