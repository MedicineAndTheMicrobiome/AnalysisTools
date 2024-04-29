#!/bin/csh

../QueryEnseml_dbSNP_by_rsID.r \
	-i bact_snps.tsv \
	-t SNP \
	-o bact_snps.wChrPos.tsv
