#!/bin/csh


../Merge_Files_By_ColumnName.pl \
	-i main.tsv \
	-k MedbioID \
	-I auxil.tsv \
	-K medbio_id \
	-o combined.tsv
