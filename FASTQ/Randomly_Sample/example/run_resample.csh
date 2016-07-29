#/bin/csh

../Randomly_Sample_from_FASTQ.pl \
	-f BadTrimming050.fastq \
	-n 10 \
	-o subsample10

../Randomly_Sample_from_FASTQ.pl \
	-f BadTrimming050.fastq \
	-n 60 \
	-o subsample60

../Randomly_Sample_from_FASTQ.pl \
	-f BadTrimming050.fastq \
	-n 25 \
	-o subsample25
