#!/bin/csh
../Fit_ALR_as_Response.r \
	-s allSamples.padded.genus.summary_table.tsv \
	-f masterfactor_file.txt \
	-M model.txt \
	-x ";"
