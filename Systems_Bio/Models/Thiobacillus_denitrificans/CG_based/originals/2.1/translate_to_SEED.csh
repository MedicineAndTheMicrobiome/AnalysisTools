#!/bin/csh

/local/devel/DAS/users/kli/SVN/DAS/SystemsBioTools/Translate_IDs/translate_reaction_formula.pl \
	-i Thiobacillus_Denitrificans_react_biggremap.tsv \
	-o Thiobacillus_Denitrificans_react.tsv \
	-m CompoundList_jiao.txt \
	-c 2 \
	-s 

/local/devel/DAS/users/kli/SVN/DAS/SystemsBioTools/Translate_IDs/translate_reaction_formula.pl \
	-i Thiobacillus_Denitrificans_react.tsv \
	-o Thiobacillus_Denitrificans_react_regress.tsv \
	-m CompoundList_jiao.txt \
	-c 2 \

