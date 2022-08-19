#!/bin/csh

../Analyze_Taxonomic_Classifications.r \
	-t read_to_group.map.unique.good.filter.unique.precluster.pick.16S_Reference.wang.taxonomy \
	#-c read_to_group.map.unique.good.filter.unique.precluster.pick.count_table\
	-c read_to_group.map.unique.good.filter.unique.precluster.pick.count_table.short \
	-o result.short
