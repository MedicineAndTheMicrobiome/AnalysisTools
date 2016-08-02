#!/bin/csh

../Estimate_SeqDepth_from_Classified_BlastHits.r \
	-b merged.blast \
	-r 10000 \
	-q 15000 \
	-o merged

