#!/bin/csh

../Calculate_Absolutes_wSpikeIn.r \
	-s example2.genus.cmF.cln.DNA.summary_table.tsv \
	-r reference_spikein.info.tsv \
	-o example.DNA

if(0)then

../Calculate_Absolutes_wSpikeIn.r \
	-s example2.genus.cmF.cln.RNA.summary_table.tsv \
	-r reference_spikein.info.tsv \
	-o example.RNA

endif
