#!/bin/csh

if (1) then
	./Find_Overexpressed_Genes.r \
		-a trans/rifle1.Thiobacillus.kegg.genes.m8blast.out.95.comp_id.summary_table.xls \
		-b metag/rifle1.Thiobacillus.kegg.genes.m8blast.out.95.comp_id.summary_table.xls \
		-o Outputfile.tbd \
		-A Transcriptomics \
		-B Genomics \
		-m rf_td_features.map.txt \
		-h dans_core/TBD_noProtSyn.core 
endif

if (1) then
	./Find_Overexpressed_Genes.r \
		-a trans/rifle1.Rhodoferax.kegg.genes.m8blast.out.95.comp_id.summary_table.xls \
		-b metag/rifle1.Rhodoferax.kegg.genes.m8blast.out.95.comp_id.summary_table.xls \
		-o Outputfile.rfr \
		-A Transcriptomics \
		-B Genomics \
		-m rf_td_features.map.txt \
		-h dans_core/RFR_noProtSyn.core 
endif
