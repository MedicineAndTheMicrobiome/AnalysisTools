#!/bin/csh

../Export_LR_Values_Batch.r \
	-t "GO.cellular_component/*names.summary_table.tsv" \
	-o GO.cellular_component
	
