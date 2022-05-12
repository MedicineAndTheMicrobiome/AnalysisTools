~/git/AnalysisTools/Read_QC_Pipeline/Run_Read_QC.pl \
	-l ./sample_list/Sample_List.FASTQ_Generation_20220308_NS_0159_Metagenome_CMM_Pool_1_2_pool12_aw_test \
	-p read_qc_30.ini \
	-r ./Raw_Data \
	-o ./QV_Data \
	-c "subsample,dust,qvtrim,seq_adapt_trim,primer_trim,contam.hs,contam.px"
