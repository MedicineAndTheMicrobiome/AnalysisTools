#!/bin/csh

../Report_Fastq_Stats_PerRead.pl \
	-f minion.o33.fastq > minion.stats.tsv

../Report_Fastq_Stats_PerRead.pl \
	-f ShortEnough.o66.fastq \
	-o 66 > ShortEnough.o66.stats.tsv
