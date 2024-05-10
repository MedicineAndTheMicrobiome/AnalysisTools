#!/bin/csh

../Extract_Variants_by_Loci.r \
	-v example_vcf.lst \
	-s targets.snps.tsv \
	-o results
