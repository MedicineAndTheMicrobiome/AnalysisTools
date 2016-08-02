#!/bin/csh

set ex = ExampleInputFiles
Map_Read_toPathway_viaEC.pl \
	-p $ex/reaction_info.tsv \
	-c $ex/pathway_type_commonname.tsv \
	-t $ex/sample_target_list.tsv \
	-o read_to_pathway
