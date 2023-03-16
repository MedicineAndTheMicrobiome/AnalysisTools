#!/bin/csh

if (0) then
	../Fit_GeneralizedNormal.r \
		-f ser_untar_metab.short100.tsv \
		-o ser_untar_metab.short100 \
		-F 

else
	../Fit_GeneralizedNormal.r \
		-f ser_untar_metab.short100.tsv \
		-o ser_untar_metab.short100 
endif
