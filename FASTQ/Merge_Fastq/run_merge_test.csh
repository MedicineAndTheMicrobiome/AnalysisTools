#!/bin/csh

Merge_Fastq.pl -f perfect/TAGGCATG-CTCTCTAT_R1.fastq -r  perfect/TAGGCATG-CTCTCTAT_R2.fastq -o perfect 
Merge_Fastq.pl -f good/TAGGCATG-CTCTCTAT_R1.75.fastq -r good/TAGGCATG-CTCTCTAT_R2.75.fastq -o good
Merge_Fastq.pl -f bad/TAGGCATG-CTCTCTAT_R1.fastq -r bad/TAGGCATG-CTCTCTAT_R2.fastq -o bad

Merge_Fastq.pl -f perfect/TAGGCATG-CTCTCTAT_R1.fastq -r  perfect/TAGGCATG-CTCTCTAT_R2.fastq -o perfect -t /usr/local/scratch/METAGENOMICS/kli
Merge_Fastq.pl -f good/TAGGCATG-CTCTCTAT_R1.75.fastq -r good/TAGGCATG-CTCTCTAT_R2.75.fastq -o good -t /usr/local/scratch/METAGENOMICS/kli
Merge_Fastq.pl -f bad/TAGGCATG-CTCTCTAT_R1.fastq -r bad/TAGGCATG-CTCTCTAT_R2.fastq -o bad -t /usr/local/scratch/METAGENOMICS/kli
