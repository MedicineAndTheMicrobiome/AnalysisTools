#!/bin/csh

../Extract_Categories_With_Table.r \
	-i Phyllo_ECIDs.cts_75.combined.main.summary_table.xls \
	-t CompatibleSolutes.tsv \
	-c 4 \
	-r 2 \
	-o extracted

../Extract_Categories_With_Table.r \
	-i Phyllo_ECIDs.cts_75.combined.main.summary_table.xls \
	-t CompatibleSolutes.tsv \
	-c 4 \
	-r 1 \
	-o extracted

../Extract_Categories_With_Table.r \
	-i Phyllo_ECIDs.cts_75.combined.main.summary_table.xls \
	-t CompatibleSolutes.tsv \
	-c 4 \
	-r 4 \
	-o extracted
