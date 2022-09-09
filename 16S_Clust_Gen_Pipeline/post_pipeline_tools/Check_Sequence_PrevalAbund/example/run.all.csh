#!/bin/csh

../Check_Sequence_PrevalAbund.r \
	-t read_to_group.map.unique.good.filter.unique.precluster.denovo.vsearch.accnos \
	-c read_to_group.map.unique.good.filter.unique.precluster.count_table \
	-o result

#read_to_group.map.unique.good.filter.unique.precluster.count_table
#read_to_group.map.unique.good.filter.unique.precluster.denovo.vsearch.accnos
#read_to_group.map.unique.good.filter.unique.precluster.fasta
