#!/bin/csh

../Build_TaxaFunct_Table.pl \
	-t taxa_list.tsv \
	-f M00165.tsv \
	-n sample_annotation.cutoff.90.tsv \
	-c 9 \
	-x sample_taxa.cutoff.90.taxa_ids.tsv \
	-o M00165



