#!/bin/csh

if(1) then
../Analyze_Fastq_QVs.r \
	-l example_trimming_list.tsv \
	-f 33 \
	-o example.wNames \
	-h 740
endif

if(1) then
../Analyze_Fastq_QVs.r \
	-l example_list.tsv \
	-f 33 \
	-o example.noNames
endif	
