#!/bin/csh

echo "Testing to 33 conversions:"
Convert_QVal_Offset.pl -f Example.64.fastq -o 64_to_33.fastq
Convert_QVal_Offset.pl -f Example.33.fastq -o 33_to_33.fastq

echo "";echo ""
echo "Testing to 64 conversions:"
Convert_QVal_Offset.pl -f Example.64.fastq -t 64 -o 64_to_64.fastq
Convert_QVal_Offset.pl -f Example.33.fastq -t 64 -o 33_to_64.fastq

echo "";echo ""
echo "Testing back conversions:"
Convert_QVal_Offset.pl -f 64_to_33.fastq -t 64 -o 64.fastq
Convert_QVal_Offset.pl -f 33_to_64.fastq -t 33 -o 33.fastq

echo "";echo ""
echo "Testing corrupted file:"
Convert_QVal_Offset.pl -f Example.mixed.fastq -t 64 -o tmp.fastq
Convert_QVal_Offset.pl -f Example.mixed.fastq -t 33 -o tmp.fastq

echo ""; echo ""
echo "Testing null conversions (Should be identical)"
diff -qs 64_to_64.fastq Example.64.fastq
diff -qs Example.33.fastq 33_to_33.fastq

echo ""
echo "Testing forward conversions (Should be different)"
diff -qs 64_to_33.fastq Example.64.fastq
diff -qs 33_to_64.fastq Example.33.fastq

echo ""
echo "Testing back conversions (Should be identical)"
diff -qs 33.fastq Example.33.fastq
diff -qs 64.fastq Example.64.fastq

rm 64_to_33.fastq 33_to_33.fastq 64_to_64.fastq 33_to_64.fastq 64.fastq 33.fastq tmp.fastq
