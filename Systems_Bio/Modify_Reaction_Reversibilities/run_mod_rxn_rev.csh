#!/bin/csh

Modify_Reaction_Reversibilities.pl \
	-i Thiobacillus_Denitrificans_react.BiGG.tsv \
	-a output.action.tsv \
	-o Thiobacillus_Denitrificans_react.BiGG.mod.tsv

Modify_Reaction_Reversibilities.pl \
	-i Thiobacillus_Denitrificans_react.tsv \
	-a output.action.tsv \
	-o Thiobacillus_Denitrificans_react.mod.tsv


