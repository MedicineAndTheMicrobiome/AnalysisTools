#!/bin/csh

# 1.) In Mothur_30/Results/Summary_Tables directory, create Collapse directory
# 2.) Copy this file into the directory
# 3.) update the st and outroot directory and run

set st=../<yourfile>.taxa.genus.cmF.cln.exp.min_0750.summary_table.tsv
set outroot=YOUR_FILENAME_ROOT

~/git/AnalysisTools/Profile/SummaryTableUtilities/Make_ReplicateTable_from_SummaryTable/Make_ReplicateTable_from_SummaryTable.r \
	-s ./$st \
	-o ./$outroot


~/git/AnalysisTools/Profile/SummaryTableUtilities/Merge_or_Subset_Samples_by_Table/Merge_or_Subset_Samples_by_Table.r \
	-i ./$st \
	-m ./$outroot\.mapping.tsv \
	-o $outroot


~/git/AnalysisTools/Column/CalculateNewColumns/calculate_new_columns.r \
	-i ./Collapsed/$outroot.Collapsed.meta.tsv \
	-f ~/git/AnalysisTools/Column/CalculateNewColumns/examples/medbio.formulas \
	-o ./Collapsed/$outroot.wSampIDForm.tsv
