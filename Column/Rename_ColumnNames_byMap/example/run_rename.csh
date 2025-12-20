#!/bin/csh

../Rename_ColumnNames_byMap.r \
	-i example.tsv \
	-m ID_to_Desc.map \
	-o output.tsv

../Rename_ColumnNames_byMap.r \
	-i example.tsv \
	-m ID_to_Desc.map \
	-o output.rfriendly.concat.tsv \
	-r -c 
