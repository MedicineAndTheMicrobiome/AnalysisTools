#!/bin/csh

../QueryEnseml_dbSNP_by_rsID.r \
	-i wang_snps.tsv \
	-t SNP \
	-o wang_snps.wChrPos.tsv
