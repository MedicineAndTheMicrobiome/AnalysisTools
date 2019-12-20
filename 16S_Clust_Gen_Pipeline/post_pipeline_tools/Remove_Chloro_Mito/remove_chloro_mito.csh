#!/bin/csh

set target_genus_summary_table=
set chloro_mito_list=~/git/AnalysisTools/16S_Clust_Gen_Pipeline/post_pipeline_tools/Remove_Chloro_Mito/chloro_mito_genus.lst

~/git/AnalysisTools/Profile/SummaryTableUtilities/Filter_Categories_By_RemoveList.r \
	-i $target_genus_summary_table \
	-l $chloro_mito_list \
	-o $target_genus_summary_table:r:r\.chlr_mito_filt.summary_table.tsv
