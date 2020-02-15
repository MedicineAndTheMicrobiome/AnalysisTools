#!/bin/csh

../Filter_Samples_by_Regex.r \
	-i canned_example.summary_table.tsv \
	-k "^00[0-9][0-9]\." \
	-o keep


../Filter_Samples_by_Regex.r \
	-i canned_example.summary_table.tsv \
	-r "^00[0-9][0-9]\." \
	-o remove

../Filter_Samples_by_Regex.r \
	-i canned_example.summary_table.tsv \
	-r "3145962654" \
	-o no_removals

../Filter_Samples_by_Regex.r \
	-i canned_example.summary_table.tsv \
	-k "2718282" \
	-o no_keepers

