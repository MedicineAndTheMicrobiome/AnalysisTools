#!/bin/csh

~/git/AnalysisTools/ProjectSpecific/Pull_16S_ForMothur/Pull_16S_ForMothur.pl \
	-d /mnt/cmmnas02/SequencingRuns \
	-s seq_run_list.lst \
	-p proj_id_list.lst \
	-o ./example_out

