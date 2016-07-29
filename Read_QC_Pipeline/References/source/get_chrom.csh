#!/bin/csh

#foreach ch (01 02 03 04 05 06 07)
#foreach ch (01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT Un)

	#@ inch = $ch
#	echo $inch $ch
	#wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_$ch/hs_ref_GRCh38.p2_chr$inch.fa.gz

#end

	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_08/hs_ref_GRCh38.p2_chr8.fa.gz
	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_09/hs_ref_GRCh38.p2_chr9.fa.gz
