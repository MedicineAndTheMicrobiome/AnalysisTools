#!/bin/csh

set DATA = /local/ifs2_scratch/METAGENOMICS/kli/subsmp/V1V3/

Extract_Sequences_For_OTU.pl \
	-o Otu00006 \
	-d ./ \
	-f /local/ifs2_scratch/METAGENOMICS/kli/subsmp/V1V3/read_to_grp.V1V3.map.fasta \
	-r $DATA/read_to_grp.V1V3.map.unique.good.filter.unique.precluster.pick.an.list \
	-p $DATA/read_to_grp.V1V3.map.unique.good.filter.unique.precluster.pick.names \
	-g $DATA/read_to_grp.V1V3.map.groups

