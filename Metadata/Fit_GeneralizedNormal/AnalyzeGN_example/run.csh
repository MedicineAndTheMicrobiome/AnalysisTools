#!/bin/csh

if(1)then

../Analyze_GenNorm_Results.r \
	-f Serum_Untargeted.transf_var.tsv \
	-b Serum_Untargeted.betas.tsv \
	-o Serum_Untargeted
else

../Analyze_GenNorm_Results.r \
	-f Saliva_Untargeted.transf_var.tsv \
	-b Saliva_Untargeted.betas.tsv \
	-o Saliva_Untargeted
endif
