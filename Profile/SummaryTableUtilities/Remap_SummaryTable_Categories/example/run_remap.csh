#!/bin/csh

../Remap_SummaryTable_Categories.r \
	-s 0159.all.cohorts.GOFuncIDs.joint.abridged.summary_table.tsv \
	-m GO_0003674.molecular_function.groupings.map \
	-n GO_0003674.molecular_function.immed_descr.tsv \
	-o GOFunct.abridged.GO_0003674
