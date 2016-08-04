#!/bin/csh

echo "Usage: $0 <sample id to filename mapping file>"

set binpath = $0:h
set rootfname = $1:t

echo "Root output name:" $rootfname

echo $binpath

$binpath/Pull_Taxonomic_Counts.pl -f $1 -p 
$binpath/Pull_Taxonomic_Counts.pl -f $1 

$binpath/Count_Unclassifiables.pl -p $rootfname\.padded
$binpath/Count_Unclassifiables.pl -p $rootfname\.unpadded

mkdir unpadded
mv $rootfname\.unpadded.* unpadded

$binpath/Graph_Taxonomic_Diversity.r -i $rootfname\.padded.domain.summary_table.xls
$binpath/Graph_Taxonomic_Diversity.r -i $rootfname\.padded.phylum.summary_table.xls
$binpath/Graph_Taxonomic_Diversity.r -i $rootfname\.padded.class.summary_table.xls
$binpath/Graph_Taxonomic_Diversity.r -i $rootfname\.padded.order.summary_table.xls
$binpath/Graph_Taxonomic_Diversity.r -i $rootfname\.padded.family.summary_table.xls
$binpath/Graph_Taxonomic_Diversity.r -i $rootfname\.padded.genus.summary_table.xls

$binpath/Rank_Abundance/RankAbundance_Analysis.r $rootfname\.padded.domain.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Analysis.r $rootfname\.padded.phylum.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Analysis.r $rootfname\.padded.class.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Analysis.r $rootfname\.padded.order.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Analysis.r $rootfname\.padded.family.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Analysis.r $rootfname\.padded.genus.summary_table.xls

$binpath/Rank_Abundance/RankAbundance_BoxPlot.r $rootfname\.padded.domain.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_BoxPlot.r $rootfname\.padded.phylum.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_BoxPlot.r $rootfname\.padded.class.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_BoxPlot.r $rootfname\.padded.order.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_BoxPlot.r $rootfname\.padded.family.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_BoxPlot.r $rootfname\.padded.genus.summary_table.xls -r

$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.domain.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.phylum.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.class.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.order.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.family.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.genus.summary_table.xls

$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.domain.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.phylum.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.class.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.order.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.family.summary_table.xls -r
$binpath/Rank_Abundance/RankAbundance_StackedBarplot.r $rootfname\.padded.genus.summary_table.xls -r

$binpath/Rank_Abundance/RankAbundance_Pie.r $rootfname\.padded.domain.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Pie.r $rootfname\.padded.phylum.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Pie.r $rootfname\.padded.class.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Pie.r $rootfname\.padded.order.summary_table.xls
$binpath/Rank_Abundance/RankAbundance_Pie.r $rootfname\.padded.family.summary_table.xls 
$binpath/Rank_Abundance/RankAbundance_Pie.r $rootfname\.padded.genus.summary_table.xls
