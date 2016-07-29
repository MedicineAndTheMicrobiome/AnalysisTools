#!/bin/csh

set reference = hs_ref_GRCh38.p2_allChr.fa
set ref_root = $reference:r:t

set tools_path = /usr/local/devel/DAS/users/kli/SVN/DAS/Read_QC_Pipeline/Dependencies

#$tools_path/bmtools/bmtools/bmtagger/bmtool -d $reference -o $ref_root.bitmask
$tools_path/bmtools/srprism mkindex -i $reference -o $ref_root.srprism
#$tools_path/blast/ncbi-blast-2.2.30+/bin/makeblastdb -in $reference -dbtype nucl
