#!/bin/csh

../../Plot_StackedBar.r \
	-i LimSC25.taxa.genus.NoControls.N500.summary_table.tsv \
	-f 0137_Lim_Metadata_20171019.tsv  \
	-o sickle.genus \
	-s ";" \
	-c "isFemale,Type"

