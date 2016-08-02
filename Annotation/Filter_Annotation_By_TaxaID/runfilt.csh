#!/bin/csh

Filter_Annotation_By_TaxaID.pl \
	-t sample/sample_taxa_summary.taxa_ids.tsv \
	-a sample/sample_annotation.cutoff.75.tsv \
	-f sample/taxa_filt.txt \
	-o sample/output_results
