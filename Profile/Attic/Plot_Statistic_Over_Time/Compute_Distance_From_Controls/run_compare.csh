#!/bin/csh

#foreach date(_0)
foreach dist(wrd euc)
	if ($dist == "euc") then
		set range = 2
	else
		set range = 5
	endif

	foreach date(_0 _1 _3 _5 _7 14)
		Compare_Two_Treatments.r \
			-t data/phaseIIFerretSampleDescription.txt \
			-c 5 \
			-d data/V1V3V3V5allSamples_samples.padded.genus.summary_table.xls.$dist.$date.distmat \
			-o data/V1V3V3V5allSamples_samples.padded.genus.summary_table.xls.$dist.$date.out \
			-r $range
	end
end
