#/bin/csh

../Extract_by_ID.pl \
	-f BadTrimming050.fastq \
	-l sample.list \
	-o extracted.fastq 

../Exclude_by_ID.pl \
	-f BadTrimming050.fastq \
	-l sample.list \
	-o exluded.fastq 

