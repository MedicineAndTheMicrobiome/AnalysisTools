#!/bin/csh

foreach level ( family order genus phylum )
./Generate_Venn_Bars.r \
	-a examples/deepSamples.$level.summary_table.xls \
	-b examples/shallowSamples.$level.summary_table.xls \
	-A Deep \
	-B Shallow \
	-o DeepVsShallow.$level 
end
