#!/bin/csh

Kmers_to_SummaryTable.pl \
	-l example/fasta_list_NOTannotated.tsv \
	-o /usr/local/scratch/METAGENOMICS/kli/kmertest \
	-k 31 \
	-p 1000
