#!/bin/csh

../MapTaxonomicAbundance_to_FBAModel.r \
	-s AQXX.kraken.bracken.joint.summary_table.tsv \
	-x g_to_model.map \
	-r override.map \
	-o AQXX
