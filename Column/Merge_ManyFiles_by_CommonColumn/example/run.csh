#!/bin/csh

../Merge_ManyFiles_by_CommonColumn.r \
	-l files.lst \
	-c MedbioID \
	-o combined.tsv
