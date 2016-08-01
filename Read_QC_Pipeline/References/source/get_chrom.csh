#!/bin/csh

# You will need to pull down the reference fasta files to perform host contaminant
# screening.  Versions keep changing so you'll need to update this script to pull
# down the latest files.

foreach ch (01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22)
	@ inch = $ch
	echo $inch $ch
	echo wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_$ch/hs_ref_GRCh38.p7_chr$inch.fa.gz
	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_$ch/hs_ref_GRCh38.p7_chr$inch.fa.gz
end

foreach ch (X Y MT Un)
	echo wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_$ch/hs_ref_GRCh38.p7_chr$ch.fa.gz
	wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/CHR_$ch/hs_ref_GRCh38.p7_chr$ch.fa.gz
end
