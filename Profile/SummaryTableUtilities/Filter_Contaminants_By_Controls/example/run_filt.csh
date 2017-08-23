#!/bin/csh

../Filter_Contaminants_By_Controls.r \
	-i grads5k.clpsd.taxa.genus.summary_table.tsv.bal_id_type.summary_table.tsv \
	-p bal_to_bctl.map \
	-l 2
