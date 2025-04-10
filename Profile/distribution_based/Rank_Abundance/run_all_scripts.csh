#!/bin/csh

echo "Usage: $0 The input file should be a summary table file as filename.>"

set binpath = $0:h
set rootfname = $1

echo "Root output name:" $rootfname

$binpath/RankAbundance_Analysis.r -i $rootfname

$binpath/RankAbundance_Analysis.r -i $rootfname -s ";"

$binpath/RankAbundance_BoxPlot.r $rootfname -r

$binpath/RankAbundance_StackedBarplot.r $rootfname

$binpath/RankAbundance_StackedBarplot.r $rootfname -r

$binpath/RankAbundance_Pie.r $rootfname
