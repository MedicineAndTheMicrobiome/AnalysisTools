#!/bin/csh

foreach dist (wrd euc)
	./Cluster_Influencers.r \
		-i grads.taxa.genus.collapsed.obsrv_cln.summary_table.tsv \
		-d $dist
end

cp grads.taxa.genus.collapsed.obsrv_cln.*.cl_inf.pdf /scratch/Temp

#../../SummaryTableUtilities/Example_SummaryTables/test5x6.summary_table.tsv
