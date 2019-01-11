#!/bin/csh

../Prescreen_Sample_woFactorNAs.r \
	-s grads.balfilt_owsh.summary_table.tsv \
	-f grads.metadata.bal.tsv \
	#-l no04S.lst \
	-S grads.stool.summary_table.tsv \
	-m grads.bal_stool.map \
	-d ./ \
	-t targeted_Var \
	-q required_var \
	-n 1000	
	
