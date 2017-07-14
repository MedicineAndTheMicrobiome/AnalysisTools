#!/bin/csh

../Generate_Summary_Tables.pl \
        -l ~/testing/Results/sample_ids.txt \
	-A ~/testing/AnnotResults \
	-S ~/testing/SummaryTables \
        -p /home/kelvinli/git/AnalysisTools/Annotation/AnnotPipeline/Assign_Annotation/annot_pipeline.ini \

