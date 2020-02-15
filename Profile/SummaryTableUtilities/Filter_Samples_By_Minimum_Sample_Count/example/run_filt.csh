#!/bin/csh

../Filter_Samples_By_Minimum_Sample_Count.r \
	-i canned_example.summary_table.tsv \
	-c 0,3000,5000,10000 \
	-p

	
