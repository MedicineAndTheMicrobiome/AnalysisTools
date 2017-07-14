#!/bin/csh

set wd = $PWD

../Generate_Summary_Tables.pl \
        -l $wd/BlastResults/sample_ids.txt \
	-A $wd/AnnotResults \
	-S $wd/SummaryTables \
	-p ~/git/AnalysisTools/Annotation/AnnotPipeline/Assign_Annotation/annot_pipeline.ini

