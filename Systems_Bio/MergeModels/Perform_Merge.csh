#!/bin/csh

set INPUT="InputModels"
set INTM="Intermediate"
set OUTP="MergedModels"

mkdir $INTM
mkdir $OUTP

dos2unix $INPUT/*

#1.) Translate external metabolite ids to BiGG, remap compartment ID to common ID

Remap_Reaction_Compounds.pl \
	-i $INPUT/Thiobacillus_Denitrificans_react.tsv \
	-o $INTM/Thiobacillus_Denitrificans_react.BiGG.tsv \
	-m CompoundList.txt \
	-c e0 \
	-n e \
	-u $INTM/Thiobacillus_Denitrificans_react.BiGG.used_ids

merge_by_map.pl \
	-i $INTM/Thiobacillus_Denitrificans_react.BiGG.used_ids \
	-m BiGG_Metabolite_Descriptions.txt \
	-k 0 \
	-s 1 \
	-o $INTM/Thiobacillus_Denitrificans_met.BiGG.tsv

#2.) Split out exchange reactions

Separate_ExchangeReactions.pl \
	-i $INTM/Thiobacillus_Denitrificans_react.BiGG.tsv \
	-o $INTM/Thiobacillus_Denitrificans_react.BiGG \

Separate_ExchangeReactions.pl \
	-i $INPUT/Rhodoferax_Ferrireducens_react.tsv \
	-o $INTM/Rhodoferax_Ferrireducens_react

#3.) Merge exchange reactions

Merge_Reactions_BasedOnStoichiometry.pl \
	-a $INTM/Rhodoferax_Ferrireducens_react.exchange.tsv \
	-b $INTM/Thiobacillus_Denitrificans_react.BiGG.exchange.tsv \
	-o $INTM/RhTh_Merged.exchange

#4.) Merge cellular reactions

foreach combo ( 20,0 19,1 18,2 17,3 16,4 15,5 14,6 13,7 12,8 11,9 10,10 9,11 8,12 7,13 6,14 5,15 4,16 3,17 2,18 1,19 0,20 )

	set A=`echo $combo | cut -f 1 -d,`
	set B=`echo $combo | cut -f 2 -d,`

	Combine_Cellular_Reactions.pl \
		-a $INTM/Rhodoferax_Ferrireducens_react.cell.tsv \
		-b $INTM/Thiobacillus_Denitrificans_react.BiGG.cell.tsv \
		-A $A \
		-B $B \
		-o $INTM/RhTh_Merged

	#5.) Concatenate merged exchanged reactions with merged cellular reactions

	cat $INTM/RhTh_Merged.$A\v$B\.tsv $INTM/RhTh_Merged.exchange > $OUTP/RhTh_$A\v$B\_react.tsv

	#6.) Concatenated metabolite file

	#cp $INPUT/Thiobacillus_Denitrificans_met.tsv $OUTP/RhTh_$A\v$B\_met.tsv
	#sed 1d $INPUT/Rhodoferax_Ferrireducens_met.tsv >> $OUTP/RhTh_$A\v$B\_met.tsv
	Combine_Metabolites.pl \
		-a $INPUT/Rhodoferax_Ferrireducens_met.tsv \
		-b $INPUT/Thiobacillus_Denitrificans_met.tsv \
		-o $INTM/RF_TD_met.BiGG.tsv

	Combine_Metabolites.pl \
		-a $INTM/RF_TD_met.BiGG.tsv \
		-b $INTM/Thiobacillus_Denitrificans_met.BiGG.tsv \
		-o $OUTP/RhTh_$A\v$B\_met.tsv

	#7.) Generate description file

	Generate_Description.pl -r $OUTP/RhTh_$A\v$B\_react.tsv -m $OUTP/RhTh_$A\v$B\_met.tsv -o $OUTP

end

