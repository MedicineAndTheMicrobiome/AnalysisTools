#!/bin/csh

if (1) then
	../Fit_GeneralizedNormal.r \
		-f ser_untar_metab.short100.tsv \
		-t ser_untar_metab.short100.targets.tsv \
		-o ser_untar_metab.short100 \
		-F 

else
	../Fit_GeneralizedNormal.r \
		-f ser_untar_metab.short100.tsv \
		-t ser_untar_metab.short100.targets.tsv \
		-o ser_untar_metab.short100 
endif
