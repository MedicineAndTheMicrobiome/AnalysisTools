#!/bin/csh

~/git/AnalysisTools/Annotation/Extract_Taxa_x_GO_Matrix/Extract_Taxa_x_GO_Matrix.r \
	-T blast.out.comp_id.trembl.taxa_est.highest.taxa_names.tsv \
	-t Genus \
	-A blast.out.comp_id.trembl.accum.40.taxa_filt.kept \
	-a GOFuncIDs \
	-m GO_0005488.binding.groupings.map \
	-o genus.binding
