#!/bin/csh

../BeforeAfterProportions.r \
	-B Before.summary_table.tsv \
	-A After.summary_table.tsv \
	-o Proportions.tsv
