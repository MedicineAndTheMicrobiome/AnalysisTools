#!/bin/csh

set list_of_targeted_variables=

cat $list_of_targeted_variables | \
	gawk '{OFS=""; print "pctd_", $1, "=perc_diff(", $1, ")"}' |
	> $list_of_targeted_variables\.formaulas
