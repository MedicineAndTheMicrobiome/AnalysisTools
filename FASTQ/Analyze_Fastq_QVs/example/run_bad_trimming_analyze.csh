#!/bin/csh

../Analyze_Fastq_QVs.r \
	-l RifleRaw_list.tsv \
	-f 33 \
	-o RifleRaw \
	-h 1000

../Analyze_Fastq_QVs.r \
	-l RifleFinalFrag_list.tsv \
	-f 33 \
	-o RifleFinalFrag \
	-h 1000

../Analyze_Fastq_QVs.r \
	-l RifleFinalPaired_list.tsv \
	-f 33 \
	-o RifleFinalPaired \
	-h 1000

../Analyze_Fastq_QVs.r \
	-l BadTrimming_list.tsv \
	-f 64 \
	-o BadTrimming
	
