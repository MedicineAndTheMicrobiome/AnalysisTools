#!/bin/csh

../Compute_Abundance_Power.r \
	-i samples.summary_table.tsv \
	-t "Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas" \
	-o Pseudomonas

../Compute_Abundance_Power.r \
	-i samples.summary_table.tsv \
	-t "Bacteria;Bacillota;Bacilli;Staphylococcales;Staphylococcaceae;Staphylococcus" \
	-o Staphylococcus



