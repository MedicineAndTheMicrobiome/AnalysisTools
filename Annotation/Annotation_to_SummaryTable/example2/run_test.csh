#!/bin/csh

mkdir results

../Annotation_to_SummaryTable.pl \
	-l file.lst \
	-o results/my_res
