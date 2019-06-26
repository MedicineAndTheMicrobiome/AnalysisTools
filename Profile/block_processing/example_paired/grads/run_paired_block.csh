#!/bin/csh

../../block_execute_models_paired.pl \
	-s grads5k.bal.summary_table.tsv \
	-S grads5k.owsh.summary_table.tsv \
	-f bal_owsh.sarc_only.metadata.tsv \
	-F bal.kit.id \
	-c covariates.txt \
	-g pfts.txt \
	-p bal_owsh.map \
	-A BAL \
	-B OWSH \
	-o ./results
