#!/bin/csh

../Clean_SummaryTable_Categories.r \
	-i gb_example.summary_table.tsv

../Clean_SummaryTable_Categories.r \
	-i canned_example.summary_table.tsv
