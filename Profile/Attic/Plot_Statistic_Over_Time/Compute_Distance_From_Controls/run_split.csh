#!/bin/csh

foreach dst (wrd euc)
	Split_Matrix_By_Treatment.r \
		-t data/phaseIIFerretSampleDescription.txt \
		-c 4 \
		-d data/V1V3V3V5allSamples_samples.padded.genus.summary_table.xls.$dst.r_distmat
end
