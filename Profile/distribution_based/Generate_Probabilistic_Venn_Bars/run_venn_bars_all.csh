#!/bin/csh

mkdir results

foreach level ( family order genus phylum )
./Generate_Probabilistic_Venn_Bars.r \
	-a examples/jamstec/deepSamples.$level.summary_table.tsv \
	-b examples/jamstec/shallowSamples.$level.summary_table.tsv \
	-A Deep \
	-B Shallow \
	-o results/DeepVsShallow.$level 
end
