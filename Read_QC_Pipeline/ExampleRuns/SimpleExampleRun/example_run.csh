#!/bin/csh

~/git/AnalysisTools/Read_QC_Pipeline/Run_Read_QC.pl \
	-l example_run.lst \
	-p example_read_qc.ini \
	-r ~/git/AnalysisTools/Read_QC_Pipeline/ExampleRuns/SimpleExampleRun/RawSequence \
	-o ~/git/AnalysisTools/Read_QC_Pipeline/ExampleRuns/SimpleExampleRun/QC \
	-c "dust,qvtrim,seq_adapt_trim,primer_trim" 


