#!/bin/csh

set pwd=`pwd`;

../../Merge_Mates.pl \
	-l fastq_list.split.tsv \
	-p \
	-r $pwd \
	-o $pwd/merge_output_split

../../Merge_Mates.pl \
	-l fastq_list.interleaved.tsv \
	-i \
	-r $pwd \
	-o $pwd/merge_output_interleaved
