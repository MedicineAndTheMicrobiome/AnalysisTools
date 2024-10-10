#!/bin/csh

foreach file (/mnt/dom_cmm_data/cmmnas02/Projects2/0159_Metagenomics/Annotation_by_type_all/function/function_all_cohorts/0159.all.cohorts.GOCompIDs.joint.summary_table.tsv /mnt/dom_cmm_data/cmmnas02/Projects2/0159_Metagenomics/Annotation_by_type_all/function/function_all_cohorts/0159.all.cohorts.GOFuncIDs.joint.summary_table.tsv /mnt/dom_cmm_data/cmmnas02/Projects2/0159_Metagenomics/Annotation_by_type_all/function/function_all_cohorts/0159.all.cohorts.GOProcIDs.joint.summary_table.tsv)
	
	echo $file:t:r:r\.abridged.summary_table.tsv
	head -n 100 $file | cut -f 1-100 > $file:t:r:r\.abridged.summary_table.tsv

end
