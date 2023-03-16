#!/bin/csh

if (0) then
	../Fit_GeneralizedNormal.r \
		-f ser_untar_metab.short101.tsv \
		-t ser_untar_metab.short101.targets.tsv \
		-o ser_untar_metab.short101 \
		-F 

else
	../Fit_GeneralizedNormal.r \
		-f ser_untar_metab.short101.tsv \
		-t ser_untar_metab.short101.targets.tsv \
		-o ser_untar_metab.short101 
endif
