#!/bin/csh

../Compute_Pathway_Likelihoods.r \
	-i imputed.ECs.abundance.summary_table.tsv \
	-o imputed.ECs.abundance \
	-p reactions.dat.ec_pathway.tsv \
	-n pathways.col


