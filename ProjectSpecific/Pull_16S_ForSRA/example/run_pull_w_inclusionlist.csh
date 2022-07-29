#!/bin/csh

~/git/AnalysisTools/ProjectSpecific/Pull_16S_ForSRA/Pull_16S_ForSRA.pl \
	-d /mnt/cmmnas02/SequencingRuns \
	-s seq_run_list.lst \
	-p proj_id_list.lst \
	-i include.lst \
	-o ./example_out_incllist \
	

