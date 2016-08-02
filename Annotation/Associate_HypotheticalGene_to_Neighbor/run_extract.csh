#!/bin/csh

Associate_HypotheticalGene_to_Neighbor.pl \
	-f example_inputs/Rhodoferax_ferrireducens_chromosome.feature_table.tsv \
	> rf_features.map.txt

Associate_HypotheticalGene_to_Neighbor.pl \
	-f example_inputs/Rhodoferax_ferrireducens_plasmid.feature_table.tsv \
	>> rf_features.map.txt

Associate_HypotheticalGene_to_Neighbor.pl \
	-f example_inputs/Thiobacillus_denitrificans_chromosome.feature_table.tsv \
	> td_features.map.txt

# Add organism prefix
sed -i 's/^Rfer_/rfr:Rfer_/' rf_features.map.txt
sed -i 's/^Tbd_/tbd:Tbd_/' td_features.map.txt
	
# Combine into single file
cat rf_features.map.txt td_features.map.txt > rf_td_features.map.txt
