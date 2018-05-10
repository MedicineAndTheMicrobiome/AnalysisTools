#!/bin/csh

set root_dir=/mnt/data/projects/0106_PACT/20180508_Analysis
set sum_tab=$root_dir/PACTLHIVCKOW_OW.taxa.genus.no_mcl.summary_table.tsv
set factors=$root_dir/pact_lhiv_carl.wform
set covar=$root_dir/covariates.txt
#set var_grp=Echocardiogram_variables.txt
set var_grp="Plasma_biomarkers.txt Echocardiogram_variables.txt CT_scan_variables.txt PulmonaryFunctionTests.txt"



foreach gr ($var_grp)

	echo working on $root_dir/$gr

	#cat $covar $root_dir/$gr > required.tmp

	~/git/AnalysisTools/Profile/distribution_based/Fit_Diversity_as_Predictor/Fit_Diversity_as_Predictor.r \
		-s $sum_tab \
		-f $factors \
		-y $root_dir/$gr \
		-o ./$gr:t\.out \
		-c $covar \
		-q $root_dir/$gr
end
	
