#!/bin/csh

set test_dir=~/git/AnalysisTools/16S_Clust_Gen_Pipeline/testing_files/

~/git/AnalysisTools/16S_Clust_Gen_Pipeline/Run_Mothur_Steps.pl \
        -f $test_dir/sample_list.fasta \
        -g $test_dir/sample_list.groups \
	-r $test_dir/abridged.reference.align \
	-o $test_dir/test_results \
	-p 16 \
	-c 0.1

