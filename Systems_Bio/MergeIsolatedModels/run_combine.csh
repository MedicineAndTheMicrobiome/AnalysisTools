#!/bin/csh

foreach cond (`ls Bounds | grep .tsv | cut -f 2 -d .`)

Combine_Models.pl \
	-a SourceModels/Rhodoferax_Ferrireducens \
	-b SourceModels/Thiobacillus_Denitrificans \
	-x Bounds/shared_exch_rxn_bounds\.$cond\.tsv \
	-A e=RFe,c=RFc\
	-B e=TDe,c=TDc\
	-e RFe \
	-f TDe \
	-i CompoundList.txt \
	-o CombinedModels/Combined\.$cond

end
