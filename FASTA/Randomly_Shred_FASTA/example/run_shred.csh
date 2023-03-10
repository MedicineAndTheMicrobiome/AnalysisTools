#!/bin/csh

../Randomly_Shred_FASTA.pl \
	-f sequence.fasta \
	-n 10 \
	-l 1000 \
	-o shredded.fasta	

../Randomly_Shred_FASTA.pl \
	-f shredded.fasta \
	-n 4 \
	-l 100 \
	-o shredded_2.fasta	

