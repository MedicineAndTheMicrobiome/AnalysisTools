#!/bin/csh

foreach level ( family order genus phylum )
./Generate_Probabilistic_Venn_Bars.r \
	-a examples/deepSamples.$level.summary_table.tsv \
	-b examples/shallowSamples.$level.summary_table.tsv \
	-A Deep \
	-B Shallow \
	-o DeepVsShallow.$level 
end
