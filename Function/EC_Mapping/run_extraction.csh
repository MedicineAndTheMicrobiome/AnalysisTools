#!/bin/csh

Extract_EC_Descriptions.pl \
	-t enzyme.dat \
	> EC_Info.tsv

Extract_EC_Class_Descriptions.pl \
	-t enzclass.txt \
	>> EC_Info.tsv
	

