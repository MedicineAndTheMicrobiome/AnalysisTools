#!/bin/csh

remap_column_values.pl -i sorted_uniref_results.btab.perc_id -m Uniprot_TaxaID.tsv -c 2 -p -f > sorted_uniref_results.btab.perc_id.taxa_id
