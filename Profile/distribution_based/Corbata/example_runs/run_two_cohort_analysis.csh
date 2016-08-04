#!/bin/csh

mkdir example_two_cohort_output

../Run_Two_SampleComparisons.pl \
	-a ../example_data/cohortA.phylum.summary_table.tsv \
	-b ../example_data/cohortB.phylum.summary_table.tsv \
	-o example_two_cohort_output/two_cohort

