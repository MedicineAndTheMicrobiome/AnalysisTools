#!/bin/csh

foreach file (`ls *.csv`)
./SNPs_to_DistMat.r \
	-s $file \
	-k target_strains.txt
end
