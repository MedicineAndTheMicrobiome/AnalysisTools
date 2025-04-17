#!/bin/csh


~/git/AnalysisTools/Profile/SummaryTableUtilities/Select_Categories_By_PCA/Select_Categories_By_PCA.r \
	-v 0.025 \
	#-s GO_0044281.small_molecule_metabolic_process.tpn100.prof.summary_table.tsv \
	#-s GO_0110165.cellular_anatomical_entity.tpn100.prof.summary_table.tsv \
	-s GO_1902494.catalytic_complex.tpn100.prof.summary_table.tsv \
	#-s GO_0032991.protein_containing_complex.tpn100.prof.summary_table.tsv \
	#-s GO_0042592.homeostatic_process.tpn100.prof.summary_table.tsv \
	-o example
