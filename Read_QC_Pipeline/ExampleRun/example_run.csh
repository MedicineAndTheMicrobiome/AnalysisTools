#!/bin/csh

../Run_Read_QC.pl \
	-l example_run.lst \
	-p example_read_qc.ini \
	-r ~/git/AnalysisTools/Read_QC_Pipeline/ExampleRun/RawSequence \
	-o ~/git/AnalysisTools/Read_QC_Pipeline/ExampleRun/QC \
	-c "contam.hs,contam.rn,contam.px,contam.ec,dust,qvtrim"


