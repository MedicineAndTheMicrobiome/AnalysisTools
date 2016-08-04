#!/bin/csh

./Compute_OverExpression.r \
	-T Transcriptomic.summary_table.tsv \
	-G Genomic.summary_table.tsv \
	-o Outputfile
