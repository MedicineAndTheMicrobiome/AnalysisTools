#!/bin/csh


foreach condition ( `ls Bounds | grep .tsv | cut -f 2 -d .` )
	Analyze_Merged_Model.r \
		-m CombinedModels/Combined\.$condition \
		-A BIO_Rfer3.RFe \
		-B ecoli_biomass.TDe \
		-a Rhodoferax_ferrireducens \
		-b Thiobacillus_denitrificans \
		-o Results/Combined\.$condition 
end
