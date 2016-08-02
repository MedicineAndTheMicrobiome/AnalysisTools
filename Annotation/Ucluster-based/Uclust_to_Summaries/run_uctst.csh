#!/bin/csh 

# Short test
if(0) then
Uclust_to_Summary_Table.r \
	-i examples/sample_list.short.tsv \
	-o examples/sample.short.out \
	-t
endif

# Complete test
if(1) then
Uclust_to_Summary_Table.r \
	-i examples/sample_list.tsv \
	-o examples/phyllo.surfactant_uclust \
	-t

endif
