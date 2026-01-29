#!/bin/csh

if(1)then

~/git/AnalysisTools/Profile/SummaryTableUtilities/Extract_Categories_By_KeywordList.r \
	-i AQXX.kraken.bracken.joint.summary_table.tsv \
	-k AQXX.minabd.0.01.map \
	-E -R \
	-o AQXX.kraken.bracken.joint.selected


~/git/AnalysisTools/Profile/SummaryTableUtilities/DoubleNormalize_noLowAbundance/DoubleNormalize_noLowAbundance.r \
	-i AQXX.kraken.bracken.joint.selected.summary_table.tsv \
	-c 0.001 \
	-o AQXX.kraken.bracken.joint.selected


~/git/AnalysisTools/Profile/SummaryTableUtilities/Rename_Summary_Table_Categories.pl \
	-i AQXX.kraken.bracken.joint.selected.mincutoff_0.001.summary_table.tsv \
	-m model_cat.map \
	-r \
	-o AQXX.kraken.bracken.joint.selected.mincutoff_0.001.xml.tsv

endif
