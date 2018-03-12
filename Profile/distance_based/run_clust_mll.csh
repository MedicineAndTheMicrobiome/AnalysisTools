#!/bin/csh


set dist="man"
set metadata=../../Metadata/20171109_OKeefe-Metadata_NOCONTROLS.tsv
set st=../../Count-Based-Analysis/SummaryTables/OKeefe25_N5000.taxa.genus.NOCONTROLS.summary_table.tsv
set out=./okeef.taxa.genus
set model="Population+Environment+Population:Environment"

~/git/AnalysisTools/Profile/distance_based/Cluster_MultinomLogLin/Cluster_MultinomLogLin.r \
	-i $st \
	-f $metadata \
	-m $model \
	-o $out \
	-d $dist 

~/git/AnalysisTools/Profile/distance_based/Plot_ClusterHeatmap/Plot_ClusterHeatmap.r \
	-i $st \
	-f $metadata \
	-d $dist \
	-m $model \
	-c $out\.$dist\.cl_mll.used_samp.lst \
	-o $out

~/git/AnalysisTools/Profile/distance_based/CH_Cluster_Stop/CH_Cluster_Stop.r \
	-i $st \
	-d $dist \
	-c $out\.$dist\.cl_mll.used_samp.lst \
	-o $out

~/git/AnalysisTools/Profile/distance_based/Cluster_Influencers/Cluster_Influencers.r \
	-i $st \
	-d $dist \
	-l $out\.$dist\.cl_mll.used_samp.lst \
	-o $out

