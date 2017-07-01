#!/bin/csh

../Generate_Alignments/Generate_Alignments.pl \
	-l ./sample_list.txt \
	-p ./test_annot_pipeline.ini \
	-o ./test_out_annotation \
	-r /home/kelvinli/git/AnalysisTools/Annotation/AnnotPipeline/Testing/fastas

../Generate_Alignments/Generate_Alignments.pl \
	-l ./sample_list.txt \
	-p ./test_frag_recruit.ini \
	-o ./test_out_fragrecruit \
	-r /home/kelvinli/git/AnalysisTools/Annotation/AnnotPipeline/Testing/fastas
