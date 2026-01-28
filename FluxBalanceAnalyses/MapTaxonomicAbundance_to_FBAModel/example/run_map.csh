#!/bin/csh

../MapTaxonomicAbundance_to_FBAModel_min_cov.r \
	-s AQXX.kraken.bracken.joint.summary_table.tsv \
	-x g_to_model.map \
	-r override.map \
	-o AQXX
