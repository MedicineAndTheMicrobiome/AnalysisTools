#!/bin/csh

set test_dir=/usr/local/devel/DAS/users/kli/SVN/DAS/16sDataAnalysis/trunk/16S_Cluster_Generation/testing_files

/usr/local/devel/DAS/users/kli/SVN/DAS/16sDataAnalysis/trunk/16S_Cluster_Generation/Run_Mothur_Steps.pl \
        -f $test_dir/sample_list.fasta \
        -g $test_dir/sample_list.groups \
	-r $test_dir/abridged.reference.align \
	-o /usr/local/scratch/METAGENOMICS/sample_list.test_dir \
	-p 16 \
	-c 0.1

