#!/bin/csh

# GO:0008150 biological_process
# GO:0016032 viral_process
# GO:0075733 intracellular transport of virus
# GO:0003674 molecular_function
# GO:0009987 cellular process

./find_GO_children.r \
	-i obo.table.tsv \
	-d GO_Groups \
	-p parent_targets.lst
