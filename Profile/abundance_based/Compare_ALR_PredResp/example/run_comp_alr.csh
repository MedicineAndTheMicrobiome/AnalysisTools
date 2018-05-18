#!/bin/csh

../Compare_ALR_PredResp.r \
	-x Plasma_biomarkers.txt.out.alr_as_pred.pvals.tsv \
	-y Plasma_biomarkers.txt.out.alr_as_resp.pvals.tsv \
	-o Plasma_biomarkers
