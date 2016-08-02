#!/bin/csh

Estimate_Taxa_from_Alignment.pl \
	-a uniref_to_percid_taxaid/sorted_uniref_results.btab.perc_id.taxa_id \
	#-a testcase \
	-t nodes.tsv \
	-n names.tsv \
	-o output.tsv

