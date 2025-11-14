#!/bin/csh

if(0)then
../reformat_olink_to_metadata.r \
	-i olink_example.raw.tsv \
	-o olink_example.bySampID.tsv
endif

../reformat_olink_to_metadata.r \
	-i olink.plasma.tsv \
	-o oline.all
