#!/bin/csh

../LassoSelect_MultinomialResponses.r \
	-f comorbcl_allsnps.tsv \
	-r response.lst \
	-c covariates.lst \
	-t Inflam27.lst \
	#-t all_snps.lst \
	-o results
