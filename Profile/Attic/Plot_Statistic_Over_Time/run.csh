#!/bin/csh

Plot_Statistic_Over_Time.r \
	-t data/phaseIIFerretSampleDescription.txt \
	-d 4 \
	-c 5 \
	-v data/V1V3V3V5allSamples_samples.padded.genus.indices.csv \
	-s 1 
	#-s 5 

Plot_Statistic_Over_Time.r \
	-t data/phaseIIFerretSampleDescription.txt \
	-d 4 \
	-c 5 \
	-v data/V1V3V3V5allSamples_samples.padded.genus.indices.csv \
	-s 5 

