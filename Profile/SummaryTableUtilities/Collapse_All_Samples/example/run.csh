#!/bin/csh

../Collapse_All_Samples.r \
	-i example.summary_table.tsv \
	-o example_default.summary_table.tsv

../Collapse_All_Samples.r \
	-i example.summary_table.tsv \
	-s new_sample_id \
	-o example_specified.summary_table.tsv

../Collapse_All_Samples.r \
	-i example.summary_table.tsv \
	-o example_filename.summary_table.tsv \
	-f 

../Collapse_All_Samples.r \
	-i example.summary_table.tsv \
	-o example_samproot.summary_table.tsv \
	-r
