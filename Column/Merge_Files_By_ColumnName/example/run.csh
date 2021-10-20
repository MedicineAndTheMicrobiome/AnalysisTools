#!/bin/csh


../Merge_Files_By_ColumnName.pl \
	-i main.tsv \
	-k MedbioID \
	-I aux.tsv \
	-K medbio_id \
	-o combined.tsv
