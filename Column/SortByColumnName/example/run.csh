#!/bin/csh

../Sort_by_ColumnName.pl \
	-i example.tsv \
	-k sample_id \
	-o example.bysample_id.tsv

../Sort_by_ColumnName.pl \
	-i example.tsv \
	-k total \
	-n \
	-o example.bytotal.tsv

../Sort_by_ColumnName.pl \
	-i example.tsv \
	-k total \
	-r \
	-n \
	-o example.bytotal_reverse.tsv

../Sort_by_ColumnName.pl \
	-i example.tsv \
	-k MedbioID \
	-o example.byMedbioID.tsv

../Sort_by_ColumnName.pl \
	-i example.tsv \
	-k total \
	-u \
	-o example.bytotalunique.tsv
