#!/bin/csh

../Plot_RankAbundance_wFactors.r \
	-i LimSC25.taxa.genus.NoControls.N500.summary_table.tsv \
	-f 0137_Lim_Metadata_20171019.tsv  \
	-T 20 \
	-o sickle.genus \
	-s ";"

../Plot_RankAbundance_wFactors.r \
	-i LimSC25.taxa.genus.NoControls.N500.summary_table.tsv \
	-f 0137_Lim_Metadata_20171019.tsv  \
	-T 10 \
	-o sickle.genus \
	-s ";"

../Plot_RankAbundance_wFactors.r \
	-i LimSC25.taxa.genus.NoControls.N500.summary_table.tsv \
	-f 0137_Lim_Metadata_20171019.tsv  \
	-T 40 \
	-o sickle.genus \
	-s ";"

