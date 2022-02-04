#!/bin/csh

../Denormalize_By_SampleID_Columns.r \
	-i example_abridge.tsv \
	-s SampleID_tumor,SampleID_normal,SampleID_stool \
	-o example_abridge.dnorm.tsv
