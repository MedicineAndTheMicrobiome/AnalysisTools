#!/bin/csh

./Compare_Taxonomic_Abundances.r \
	-a rifle1.summary_table.trans.csv \
	-b rifle1.summary_table.genom.csv \
	-o Outputfile \
	-A Transcriptomics \
	-B Genomics 
