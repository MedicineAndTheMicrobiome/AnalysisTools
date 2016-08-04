#!/bin/csh

mkdir example_single_cohort_output

../Run_Single_SampleAnalyses.pl \
	-s ../example_data/cohortA.phylum.summary_table.tsv \
	-o example_single_cohort_output/cohortA

