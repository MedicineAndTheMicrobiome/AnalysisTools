#!/bin/csh


set ref_root=hs_ref_all
set tools_path=../Dependencies

#
# 0. First download reference.  Run source/get_chrom.csh
#

cd source
get_chrom.csh
cd ..

#
# 1. Concatenate references
#

zcat source/hs_ref_GRCh38.p7_chr*fa.gz > $ref_root.fa

#
# 2. Create bitmask file
#

$tools_path/bmtools/bmtools/bmtagger/bmtool -d $ref_root.fa -o $ref_root.bitmask

#
# 3. Create srprism file
#

$tools_path/bmtools/srprism mkindex -i $ref_root.fa -o $ref_root.srprism

#
# 4. Create blast db file
#

$tools_path/blast/ncbi-blast-2.2.30+/bin/makeblastdb -in $ref_root.fa -dbtype nucl



