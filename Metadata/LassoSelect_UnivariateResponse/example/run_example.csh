#!/bin/csh

../LassoSelect_UnivariateResponse.r \
	-f Telomere_Dev_SNPs_comorb.tsv \
	-r telom_obslen_mean \
	-c covariates.lst \
	-t Telomere_SNPs.lst \
	-o results
