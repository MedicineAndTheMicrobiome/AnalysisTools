#!/bin/csh

set binpath = $0:h

foreach level (phylum class order family genus)
foreach level(main)
	
	set target=$1.padded.$level.summary_table.xls

	$binpath/RankAbundance_Analysis.r $target
	$binpath/RankAbundance_BoxPlot.r $target
	$binpath/RankAbundance_Pie.r $target
	$binpath/RankAbundance_StackedBarplot.r $target

end


