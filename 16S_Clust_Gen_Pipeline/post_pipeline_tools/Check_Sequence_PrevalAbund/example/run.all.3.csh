#!/bin/csh

../Check_Sequence_PrevalAbund.r \
	-t 0102_ALIR_2020.unique.good.filter.unique.precluster.denovo.vsearch.accnos \
	-c 0102_ALIR_2020.unique.good.filter.unique.precluster.count_table \
	-f /mnt/cmmnas02/Projects2/0102_ALIR_2020/Mothur_25/Results/0102_ALIR_2020.unique.good.filter.unique.precluster.fasta \
	-o 0102_ALIR_2020b

#read_to_group.map.unique.good.filter.unique.precluster.count_table
#read_to_group.map.unique.good.filter.unique.precluster.denovo.vsearch.accnos
#read_to_group.map.unique.good.filter.unique.precluster.fasta
