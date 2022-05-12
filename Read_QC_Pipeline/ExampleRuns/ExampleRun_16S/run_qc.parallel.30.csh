#!/bin/csh

set multiplier=4
mkdir QC

foreach offset ( `seq 1 $multiplier` )

~/git/AnalysisTools/Read_QC_Pipeline/Run_Read_QC.pl \
	-l Sample_Fastq_List.txt \
	-p read_qc_30.ini \
	-r ~/git/AnalysisTools/Read_QC_Pipeline/ExampleRuns/ExampleRun_16S \
	-o QC/Results \
	-c "dust,qvtrim,seq_adapt_trim,primer_trim" \
	-s $offset,$multiplier &

end
