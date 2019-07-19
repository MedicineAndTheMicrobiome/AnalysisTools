#!/bin/csh

../Plot_MDS_wFactors.r \
	-i LimSC25.taxa.genus.NoControls.N500.summary_table.tsv \
	-f 0137_Lim_Metadata_20171019.tsv  \
	-d man \
	-M model \
	-o sickle.genus 
if(0)then

../Plot_MDS_wFactors.r \
	-i LimSC25.taxa.genus.NoControls.N500.summary_table.tsv \
	-f 0137_Lim_Metadata_20171019.tsv  \
	-o sickle.genus 

../Plot_MDS_wFactors.r \
	-i LimSC25.taxa.genus.NoControls.N500.summary_table.tsv \
	-f 0137_Lim_Metadata_20171019.tsv  \
	-o sickle.genus 

endif
