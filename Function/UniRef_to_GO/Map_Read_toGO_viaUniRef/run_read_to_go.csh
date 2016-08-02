#!/bin/csh

set ex = ExampleInputFiles

Map_Read_toGO_viaUniRef.pl \
	-r $ex/reads_to_uniref.tsv.one \
	-u $ex/uniref_to_go.map.short \
	-g $ex/GO.map \
	-o read_to_go_path.test



Map_Read_toGO_viaUniRef.pl \
	-r $ex/reads_to_uniref.tsv \
	-u $ex/uniref_to_go.map.full \
	-g $ex/GO.map \
	-o read_to_go_path.example
