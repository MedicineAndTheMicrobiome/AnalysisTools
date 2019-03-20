#!/bin/csh

set tmp_fname=chloro_mito.lst
set st_name=<>.taxa.genus.summary_table.tsv

echo "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria;Mitochondria_unclassified" > $tmp_fname
echo "Bacteria;Cyanobacteria;Chloroplast;Chloroplast_unclassified;Chloroplast_unclassified;Chloroplast_unclassified" >> $tmp_fname

~/git/AnalysisTools/Profile/SummaryTableUtilities/Filter_Categories_By_RemoveList.r \
	-i $st_name \
	-l $tmp_fname \
	-o $st_name:r:r.no_chlr_mito.summary_table.tsv

