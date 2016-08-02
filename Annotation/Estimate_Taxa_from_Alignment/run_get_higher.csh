#!/bin/csh

Get_Higher_Levels.pl \
	-t output.best_taxa_id.tsv \
	-m names.tsv \
	-d nodes.tsv \
	-l levels.tsv \
	-o output.upper_levels
