#!/bin/csh

set wd = $PWD

../Run_Read_Annotation.pl \
	-l $wd/BlastResults/sample_ids.txt  \
	-P annot_test \
	-I $wd/BlastResults \
	-O $wd/AnnotResults \
	-p ~/git/AnalysisTools/Annotation/AnnotPipeline/Assign_Annotation/annot_pipeline.ini \
	-c $wd/BlastResults/contam_id.txt
	
